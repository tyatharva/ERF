#include <ERF_MM5.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>

using namespace amrex;

/* Initialize lsm data structures */
void
MM5::Init (const int& /*lev*/,
           const MultiFab& cons_in,
           const Geometry& geom,
           const Real& dt)
{
    m_dt = dt;
    m_geom = geom;

    // erf.mm5.soil_theta sets BOTH the uniform column init AND the
    // permanent Dirichlet anchor at the column bottom (the bottom ghost
    // is set once here and never updated by the diffusion update). The
    // 300 K default is 12-15 K warm for a SoCal January: with the
    // 3-m column's fundamental-mode timescale ~4H^2/(pi^2 D) ~ 266 days,
    // a wrong anchor is a year-long soil heat source, not a spin-up
    // transient. erf.mm5.diag_interval > 0 prints soil-theta stats at
    // several depths every N LSM advances (spin-up measurement).
    {
        ParmParse pp("erf.mm5");
        pp.query("soil_theta", m_theta_dir);
        pp.query("diag_interval", m_diag_interval);
    }

    Box domain = geom.Domain();
    khi_lsm    = domain.smallEnd(2) - 1;

    LsmVarMap.resize(m_lsm_size);
    LsmVarMap = {LsmVar_MM5::theta};

    LsmVarName.resize(m_lsm_size);
    LsmVarName = {"theta"};

    // NOTE: All boxes in ba extend from zlo to zhi, so this transform is valid.
    //       If that were to change, the dm and new ba are no longer valid and
    //       direct copying between lsm data/flux vars cannot be done in a parfor.

    // Set box array for lsm data
    IntVect ng(0,0,1);
    BoxArray ba = cons_in.boxArray();
    DistributionMapping dm = cons_in.DistributionMap();
    BoxList bl_lsm = ba.boxList();
    for (auto& b : bl_lsm) {
        b.setBig(2,khi_lsm);                  // First point below the surface
        b.setSmall(2,khi_lsm - m_nz_lsm + 1); // Last point below the surface
    }
    BoxArray ba_lsm(std::move(bl_lsm));

    // Set up lsm geometry
    const RealBox& dom_rb = m_geom.ProbDomain();
    const Real*    dom_dx = m_geom.CellSize();
    RealBox lsm_rb = dom_rb;
    Real lsm_dx[AMREX_SPACEDIM] = {AMREX_D_DECL(dom_dx[0],dom_dx[1],m_dz_lsm)};
    Real lsm_z_hi = dom_rb.lo(2);
    Real lsm_z_lo = lsm_z_hi - Real(m_nz_lsm)*lsm_dx[2];
    lsm_rb.setHi(2,lsm_z_hi); lsm_rb.setLo(2,lsm_z_lo);
    m_lsm_geom.define( ba_lsm.minimalBox(), lsm_rb, m_geom.Coord(), m_geom.isPeriodic() );

    // Create the data and fluxes
    for (auto ivar = 0; ivar < LsmVar_MM5::NumVars; ++ivar) {
        // State vars are CC
        Real theta_0 = m_theta_dir;
        lsm_fab_vars[ivar] = std::make_shared<MultiFab>(ba_lsm, dm, 1, ng);
        lsm_fab_vars[ivar]->setVal(theta_0);

        // Fluxes are nodal in z
        lsm_fab_flux[ivar] = std::make_shared<MultiFab>(convert(ba_lsm, IntVect(0,0,1)), dm, 1, IntVect(0,0,0));
        lsm_fab_flux[ivar]->setVal(0.);
    }
}

/* Extrapolate surface temperature and store in ghost cell */
void
MM5::ComputeTsurf ()
{
    // Expose for GPU copy
    int khi = khi_lsm;

    for ( MFIter mfi(*(lsm_fab_vars[LsmVar_MM5::theta])); mfi.isValid(); ++mfi) {
        auto box2d = mfi.tilebox(); box2d.makeSlab(2,khi);

        auto theta_array = lsm_fab_vars[LsmVar_MM5::theta]->array(mfi);

        ParallelFor( box2d, [=] AMREX_GPU_DEVICE (int i, int j, int )
        {
            theta_array(i,j,khi+1) = Real(1.5)*theta_array(i,j,khi) - myhalf*theta_array(i,j,khi-1);
        });
    }
}

/* Compute the diffusive fluxes */
void
MM5::ComputeFluxes ()
{
    // Expose for GPU copy
    int khi = khi_lsm;
    Real Dsoil = m_d_soil;
    Real dzInv = m_lsm_geom.InvCellSize(2);

    for ( MFIter mfi(*(lsm_fab_flux[LsmVar_MM5::theta])); mfi.isValid(); ++mfi) {
        auto box3d = mfi.tilebox();

        // Do not overwrite the flux at the top (comes from MOST BC)
        if (box3d.bigEnd(2) == khi+1) box3d.setBig(2,khi);

        auto theta_array = lsm_fab_vars[LsmVar_MM5::theta]->array(mfi);
        auto theta_flux  = lsm_fab_flux[LsmVar_MM5::theta]->array(mfi);

        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            theta_flux(i,j,k) = Dsoil * ( theta_array(i,j,k) - theta_array(i,j,k-1) ) * dzInv;
        });
    }
}

/* Advance the solution with a simple explicit update (should use tridiagonal solve) */
void
MM5::AdvanceMM5 ()
{
    // Spin-up diagnostic: soil-theta min/mean/max at several depths.
    // NOTE ocean columns receive no MOST top flux and sit at the init
    // value forever -- they dilute the mean; min/max track the land
    // signal.
    if (m_diag_interval > 0 && (m_diag_count++ % m_diag_interval == 0)) {
        MultiFab& th = *(lsm_fab_vars[LsmVar_MM5::theta]);
        for (int koff : {0, 4, 14, 29}) {
            const int k = khi_lsm - koff;
            Real mn = std::numeric_limits<Real>::max();
            Real mx = std::numeric_limits<Real>::lowest();
            Real sm = 0.0; long np = 0;
            for (MFIter mfi(th); mfi.isValid(); ++mfi) {
                Box b = mfi.validbox(); b.setRange(2, k, 1);
                if (!b.ok()) continue;
                auto const& a = th.const_array(mfi);
                Real bmn, bmx, bsm;
                {
                    ReduceOps<ReduceOpMin, ReduceOpMax, ReduceOpSum> rops;
                    ReduceData<Real, Real, Real> rdata(rops);
                    rops.eval(b, rdata, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
                              -> GpuTuple<Real, Real, Real>
                    { return {a(i,j,kk), a(i,j,kk), a(i,j,kk)}; });
                    auto tup = rdata.value(rops);
                    bmn = get<0>(tup); bmx = get<1>(tup); bsm = get<2>(tup);
                }
                mn = std::min(mn, bmn); mx = std::max(mx, bmx);
                sm += bsm; np += b.numPts();
            }
            ParallelDescriptor::ReduceRealMin(mn);
            ParallelDescriptor::ReduceRealMax(mx);
            ParallelDescriptor::ReduceRealSum(sm);
            ParallelDescriptor::ReduceLongSum(np);
            Print() << "MM5SOIL depth " << (koff + 1) * m_dz_lsm
                    << " m: min " << mn << " mean " << sm / np
                    << " max " << mx << std::endl;
        }
    }
    // Expose for GPU copy
    Real dt = m_dt;
    Real dzInv = m_lsm_geom.InvCellSize(2);

    for ( MFIter mfi(*(lsm_fab_vars[LsmVar_MM5::theta])); mfi.isValid(); ++mfi) {
        auto box3d = mfi.tilebox();

        auto theta_array = lsm_fab_vars[LsmVar_MM5::theta]->array(mfi);
        auto theta_flux  = lsm_fab_flux[LsmVar_MM5::theta]->array(mfi);

        ParallelFor( box3d, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            theta_array(i,j,k) += dt * ( theta_flux(i,j,k+1) - theta_flux(i,j,k) ) * dzInv;
        });
    }
}
