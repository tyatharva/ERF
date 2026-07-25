#include "ERF_Utils.H"

using namespace amrex;

PhysBCFunctNoOp void_bc;

/**
 * Get the boxes for looping over interior/exterior ghost cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post.
 *
 * @param[in] bx box to intersect with 4 halo regions
 * @param[in] domain box of the whole domain
 * @param[in] width number of cells in (relaxation+specified) zone
 * @param[in] set_width number of cells in (specified) zone
 * @param[out] bx_xlo halo box at x_lo boundary
 * @param[out] bx_xhi halo box at x_hi boundary
 * @param[out] bx_ylo halo box at y_lo boundary
 * @param[out] bx_yhi halo box at y_hi boundary
 * @param[in] ng_vect number of ghost cells in each direction
 * @param[in] get_int_ng flag to get ghost cells inside the domain
 */
void
realbdy_interior_bxs_xy (const Box& bx,
                         const Box& domain,
                         const int& width,
                         Box& bx_xlo,
                         Box& bx_xhi,
                         Box& bx_ylo,
                         Box& bx_yhi,
                         const IntVect& ng_vect,
                         const bool get_int_ng)
{
    AMREX_ALWAYS_ASSERT(bx.ixType() == domain.ixType());

    //==================================================================
    // NOTE: X-face boxes take ownership of the overlapping region.
    //       With exterior ghost cells (ng_vect != 0), the x-face
    //       boxes will have exterior ghost cells in both x & y.
    //==================================================================

    // Domain bounds without ghost cells
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // Four boxes matching the domain
    Box gdom_xlo(domain); Box gdom_xhi(domain);
    Box gdom_ylo(domain); Box gdom_yhi(domain);

    // Trim the boxes to only include internal ghost cells
    gdom_xlo.setBig(0,dom_lo.x+width-1); gdom_xhi.setSmall(0,dom_hi.x-width+1);
    gdom_ylo.setBig(1,dom_lo.y+width-1); gdom_yhi.setSmall(1,dom_hi.y-width+1);

    // Remove overlapping corners from y-face boxes
    gdom_ylo.setSmall(0,gdom_xlo.bigEnd(0)+1); gdom_ylo.setBig(0,gdom_xhi.smallEnd(0)-1);
    gdom_yhi.setSmall(0,gdom_xlo.bigEnd(0)+1); gdom_yhi.setBig(0,gdom_xhi.smallEnd(0)-1);

    // Grow boxes to get external ghost cells only
    gdom_xlo.growLo(0,ng_vect[0]); gdom_xhi.growHi(0,ng_vect[0]);
    gdom_xlo.grow  (1,ng_vect[1]); gdom_xhi.grow  (1,ng_vect[1]);
    gdom_ylo.growLo(1,ng_vect[1]); gdom_yhi.growHi(1,ng_vect[1]);

    // Grow boxes to get internal ghost cells
    if (get_int_ng) {
        gdom_xlo.growHi(0,ng_vect[0]); gdom_xhi.growLo(0,ng_vect[0]);
        gdom_ylo.grow  (0,ng_vect[0]); gdom_yhi.grow  (0,ng_vect[0]);
        gdom_ylo.growHi(1,ng_vect[1]); gdom_yhi.growLo(1,ng_vect[1]);
    }

    // Populate everything
    bx_xlo = (bx & gdom_xlo);
    bx_xhi = (bx & gdom_xhi);
    bx_ylo = (bx & gdom_ylo);
    bx_yhi = (bx & gdom_yhi);
}


/**
 * Get the boxes for looping over interior/exterior ghost cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post.
 *
 * @param[in] bx box to intersect with 4 halo regions
 * @param[in] domain box of the whole domain
 * @param[in] width number of cells in (relaxation+specified) zone
 * @param[in] set_width number of cells in (specified) zone
 * @param[out] bx_xlo halo box at x_lo boundary
 * @param[out] bx_xhi halo box at x_hi boundary
 * @param[out] bx_ylo halo box at y_lo boundary
 * @param[out] bx_yhi halo box at y_hi boundary
 * @param[in] ng_vect number of ghost cells in each direction
 * @param[in] get_int_ng flag to get ghost cells inside the domain
 */
void
realbdy_bc_bxs_xy (const Box& bx,
                   const Box& domain,
                   const int& set_width,
                   Box& bx_xlo,
                   Box& bx_xhi,
                   Box& bx_ylo,
                   Box& bx_yhi,
                   const IntVect& ng_vect)
{
    AMREX_ALWAYS_ASSERT(bx.ixType() == domain.ixType());

    // Domain bounds without ghost cells
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    // Four boxes matching the domain
    Box gdom_xlo(domain); Box gdom_xhi(domain);
    Box gdom_ylo(domain); Box gdom_yhi(domain);

    // Get offsets from box index type
    IntVect iv_type = bx.ixType().toIntVect();
    int offx = (iv_type[0]==1) ? 0 : -1;
    int offy = (iv_type[1]==1) ? 0 : -1;

    // Stagger the boxes based upon index type
    gdom_xlo += IntVect(offx,0,0); gdom_xhi += IntVect(-offx,0,0);
    gdom_ylo += IntVect(0,offy,0); gdom_yhi += IntVect(0,-offy,0);

    // Trim the boxes to only include internal ghost cells
    gdom_xlo.setBig(0,dom_lo.x+set_width+offx-1); gdom_xhi.setSmall(0,dom_hi.x-set_width-offx+1);
    gdom_ylo.setBig(1,dom_lo.y+set_width+offy-1); gdom_yhi.setSmall(1,dom_hi.y-set_width-offy+1);

    // Remove overlapping corners from y-face boxes
    gdom_ylo.setSmall(0,gdom_xlo.bigEnd(0)+1); gdom_ylo.setBig(0,gdom_xhi.smallEnd(0)-1);
    gdom_yhi.setSmall(0,gdom_xlo.bigEnd(0)+1); gdom_yhi.setBig(0,gdom_xhi.smallEnd(0)-1);

    // Grow boxes to get external ghost cells only
    gdom_xlo.growLo(0,ng_vect[0]+offx); gdom_xhi.growHi(0,ng_vect[0]+offx);
    gdom_xlo.grow  (1,ng_vect[1]     ); gdom_xhi.grow  (1,ng_vect[1]     );
    gdom_ylo.growLo(1,ng_vect[1]+offy); gdom_yhi.growHi(1,ng_vect[1]+offy);

    // Populate everything
    bx_xlo = (bx & gdom_xlo);
    bx_xhi = (bx & gdom_xhi);
    bx_ylo = (bx & gdom_ylo);
    bx_yhi = (bx & gdom_yhi);
}


/**
 * Compute the RHS in the relaxation zone
 *
 * @param[in] time              current (total) time
 * @param[in] delta_t           timestep
 * @param[in] start_bdy_time    full time of the first time slice of boundary data
 * @param[in] final_bdy_time    full time of the  last time slice of boundary data
 * @param[in] bdy_time_interval time interval between boundary condition time stamps
 * @param[in] width             number of cells in (relaxation+specified) zone
 * @param[in] set_width         number of cells in (specified) zone
 * @param[in] geom              container for geometric information
 * @param[out] S_rhs            RHS to be computed here
 * @param[in] S_data            current value of the solution
 * @param[in] bdy_data_xlo boundary data on interior of low x-face
 * @param[in] bdy_data_xhi boundary data on interior of high x-face
 * @param[in] bdy_data_ylo boundary data on interior of low y-face
 * @param[in] bdy_data_yhi boundary data on interior of high y-face
 */
void
compute_wall_flux_correction (const Geometry& geom,
                              const MultiFab& cons_model,
                              const MultiFab& cons_tgt,
                              Real tau,
                              MultiFab& dvel)
{
    BL_PROFILE("compute_wall_flux_correction()");

    dvel.setVal(0.0);

    const Box& domain = geom.Domain();
    const auto dom_lo = lbound(domain);
    const auto dom_hi = ubound(domain);
    const Real dx     = geom.CellSize(0);
    const Real dy     = geom.CellSize(1);
    const Real dz     = geom.CellSize(2);

    for (MFIter mfi(dvel); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        // The column sums require the whole column on one box.
        if (bx.smallEnd(2) != dom_lo.z || bx.bigEnd(2) != dom_hi.z) { continue; }

        const Box bx2 = makeSlab(bx,2,bx.smallEnd(2));
        const Array4<const Real>& rm = cons_model.const_array(mfi);
        const Array4<const Real>& rt = cons_tgt.const_array(mfi);
        const Array4<Real>&       dv = dvel.array(mfi);

        const int klo = dom_lo.z, khi = dom_hi.z;

        ParallelFor(bx2, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            // Only columns that own a domain-wall face are corrected.
            const bool wx = (i == dom_lo.x) || (i == dom_hi.x);
            const bool wy = (j == dom_lo.y) || (j == dom_hi.y);
            if (!wx && !wy) { return; }

            Real M = Real(0.0), Mt = Real(0.0);
            for (int k = klo; k <= khi; ++k) {
                M  += rm(i,j,k,Rho_comp) * dz;
                Mt += rt(i,j,k,Rho_comp) * dz;
            }
            if (M <= Real(0.0)) { return; }

            // Corner columns own an x-wall and a y-wall; splitting the
            // correction between them avoids double counting.
            const Real split = (wx && wy) ? Real(0.5) : Real(1.0);
            const Real len   = wx ? dx : dy;
            dv(i,j,klo) = split * len * (Mt/M - Real(1.0)) / tau;
        });

        // Barotropic: broadcast the k=klo value up the column so consumers can
        // index it at any k.
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dv(i,j,k) = dv(i,j,klo);
        });
    }
    dvel.FillBoundary(geom.periodicity());
}

void
realbdy_compute_interior_ghost_rhs (const Real& time,
                                    const Real& delta_t,
                                    const Real& start_bdy_time,
                                    const Real& final_bdy_time,
                                    const Real& bdy_time_interval,
                                    const Real& nudge_factor,
                                    int width,
                                    const Geometry& geom,
                                    Vector<MultiFab>& S_rhs,
                                    Vector<MultiFab>& S_cur_data,
                                    const MultiFab* fcons_tgt,
                                    Vector<Vector<FArrayBox>>& bdy_data_xlo,
                                    Vector<Vector<FArrayBox>>& bdy_data_xhi,
                                    Vector<Vector<FArrayBox>>& bdy_data_ylo,
                                    Vector<Vector<FArrayBox>>& bdy_data_yhi,
                                    std::unique_ptr<ReadBndryPlanes>& m_r2d)
{
    BL_PROFILE_REGION("realbdy_compute_interior_ghost_RHS()");

    //
    // Note that time (= start_time+old_stage_time)  is measured as total time
    //           start_bdy_time and final_bdy_time are also measured as total time
    //

    // Get bndry data if we have it
    Vector<int>  bnd_map = {BCVars::xvel_bc, BCVars::yvel_bc, BCVars::RhoTheta_bc_comp,
                            BCVars::RhoQ1_bc_comp};
    Array4<Real> bdatxlo, bdatxhi, bdatylo, bdatyhi;
    Array4<Real> btenxlo, btenxhi, btenylo, btenyhi;
    if (m_r2d) {
        // Index is [plane orientation] and [level]
        Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time);
        bdatxlo = (*bndry_data[0])[0].array();
        bdatylo = (*bndry_data[1])[0].array();
        bdatxhi = (*bndry_data[3])[0].array();
        bdatyhi = (*bndry_data[4])[0].array();

        Vector<std::unique_ptr<PlaneVector>>& bndry_tend = m_r2d->get_tendency(time);
        btenxlo = (*bndry_tend[0])[0].array();
        btenylo = (*bndry_tend[1])[0].array();
        btenxhi = (*bndry_tend[3])[0].array();
        btenyhi = (*bndry_tend[4])[0].array();
    }

    // Relaxation constants
    Real F1 = one/(nudge_factor*delta_t);

    // Time interpolation
    Real dT = bdy_time_interval;

    int n_time    = static_cast<int>( (time-start_bdy_time) /  dT);
    int n_time_p1 = n_time + 1;
    Real alpha    = ((time-start_bdy_time) - n_time * dT) / dT;

    // Do not over run the last bdy file
    if (time >= final_bdy_time) {
        n_time    = static_cast<int>( (final_bdy_time - start_bdy_time)/ dT);
        n_time_p1 = n_time;
        alpha     = zero;
    }

    AMREX_ALWAYS_ASSERT( alpha >= zero && alpha <= one);
    Real oma   = one - alpha;

    // Temporary FABs for storage (owned/filled on all ranks)
    FArrayBox U_xlo, U_xhi, U_ylo, U_yhi;
    FArrayBox V_xlo, V_xhi, V_ylo, V_yhi;
    FArrayBox T_xlo, T_xhi, T_ylo, T_yhi;
    FArrayBox Q_xlo, Q_xhi, Q_ylo, Q_yhi;
    FArrayBox R_xlo, R_xhi, R_ylo, R_yhi;

    // Variable index map (WRFBdyVars -> Vars)
    Vector<int> var_map  = {Vars::xvel,    Vars::yvel,    Vars::cons,    Vars::cons,    Vars::cons   };
    Vector<int> ivar_map = {IntVars::xmom, IntVars::ymom, IntVars::cons, IntVars::cons, IntVars::cons};

    // Variable icomp map
    Vector<int> comp_map = {0, 0, RhoTheta_comp, RhoQ1_comp, Rho_comp};

    // Indices
    int  ivarU  = RealBdyVars::U;
    int  ivarV  = RealBdyVars::V;
    int  ivarT  = RealBdyVars::T;
    int  ivarQV = RealBdyVars::QV;
    int  ivarR  = RealBdyVars::NumTypes;   // = HindcastBdyVars::RHO; not a bdy plane

    // MASS-CONSISTENT LATERAL FORCING (erf.hindcast_mass_consistent_bdy).
    //
    // Omega is zero at k = klo (hard-zeroed) and at the SlipWall lid, so the
    // vertical mass flux telescopes out of a column sum of the discrete
    // continuity equation: the column mass tendency is determined ENTIRELY by
    // the horizontal fluxes. Prescribing rho*u and rho*v in the band therefore
    // already prescribes the column mass tendency, while rho itself is relaxed
    // toward nothing -- so the level-by-level residual has only Delta_z(rho
    // Omega) to go into, and Omega (hence w) becomes the residual variable.
    //
    // The fix is to make the target a COMPLETE state that solves ERF's own
    // discrete continuity equation:
    //   - relax Rho_comp toward the ERA5 density rho* (it is otherwise the
    //     only prognostic field in the band with no constraint at all);
    //   - build the momentum target from rho* u* rather than rho_model u*,
    //     which removes the spurious u* Delta_x(rho_model) term -- the model's
    //     own density structure entering the target's divergence as a mass
    //     source.
    // w is deliberately NOT specified: Omega comes out right as a consequence
    // of the mass budget closing, not by imposition.
    static const bool mass_consistent = [] {
        bool b = false; amrex::ParmParse pp("erf");
        pp.query("hindcast_mass_consistent_bdy", b); return b;
    }();
    const bool l_mass_consistent = mass_consistent && (fcons_tgt != nullptr);
    // Relax QV too when the state carries moisture and the boundary planes
    // include it (WRF relaxes moisture in the zone; previously only U/V/T
    // were nudged while QV was set in the specified cells only, leaving a
    // sharp qv step at the specified/relaxation interface).
    bool l_relax_qv = (S_cur_data[IntVars::cons].nComp() > RhoQ1_comp) &&
                      (!bdy_data_xlo.empty()) &&
                      (static_cast<int>(bdy_data_xlo[0].size()) > RealBdyVars::QV);
    int BdyEnd = l_relax_qv ? RealBdyVars::NumTypes : RealBdyVars::NumTypes-1;
    // Rho is appended past the bdy-plane vars; its target comes from
    // fcons_tgt, not from a boundary plane. QV is skipped inside the loops
    // when the state has none.
    //
    // erf.hindcast_relax_rho separates the two halves of the mass-consistent
    // forcing for A/B: false keeps the rho* u* momentum target but drops the
    // rho relaxation itself.
    static const bool relax_rho_bdy = [] {
        bool b = true; amrex::ParmParse pp("erf");
        pp.query("hindcast_relax_rho", b); return b;
    }();
    if (l_mass_consistent && relax_rho_bdy) { BdyEnd = ivarR + 1; }


    // NOTE: The sizing of the temporary BDY FABS is
    //       GLOBAL and occurs over the entire BDY region.

    // Size the FABs
    //==========================================================
    for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
        if (ivar == ivarQV && !l_relax_qv) { continue; }
        int ivar_idx = var_map[ivar];
        Box domain   = geom.Domain();
        auto ixtype  = S_cur_data[ivar_idx].boxArray().ixType();
        domain.convert(ixtype);

        // NOTE: Ghost cells needed for idx type mismatch between mask and data (do_upwind)
        IntVect ng_vect(0);
        //IntVect ng_vect(1,1,0);
        Box gdom(domain); gdom.grow(ng_vect);
        Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
        realbdy_interior_bxs_xy(gdom, domain, width,
                                bx_xlo, bx_xhi,
                                bx_ylo, bx_yhi,
                                ng_vect, true);

        // Size the FABs
        if (ivar  == ivarU) {
            U_xlo.resize(bx_xlo,1,The_Async_Arena()); U_xhi.resize(bx_xhi,1,The_Async_Arena());
            U_ylo.resize(bx_ylo,1,The_Async_Arena()); U_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarV) {
            V_xlo.resize(bx_xlo,1,The_Async_Arena()); V_xhi.resize(bx_xhi,1,The_Async_Arena());
            V_ylo.resize(bx_ylo,1,The_Async_Arena()); V_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarT){
            T_xlo.resize(bx_xlo,1,The_Async_Arena()); T_xhi.resize(bx_xhi,1,The_Async_Arena());
            T_ylo.resize(bx_ylo,1,The_Async_Arena()); T_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarQV){
            Q_xlo.resize(bx_xlo,1,The_Async_Arena()); Q_xhi.resize(bx_xhi,1,The_Async_Arena());
            Q_ylo.resize(bx_ylo,1,The_Async_Arena()); Q_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else if (ivar  == ivarR){
            R_xlo.resize(bx_xlo,1,The_Async_Arena()); R_xhi.resize(bx_xhi,1,The_Async_Arena());
            R_ylo.resize(bx_ylo,1,The_Async_Arena()); R_yhi.resize(bx_yhi,1,The_Async_Arena());
        } else {
            continue;
        }
    } // ivar


    // NOTE: These operations use the BDY FABS and RHO. The
    //       use of RHO to go from PRIM -> CONS requires that
    //       these operations be LOCAL. So we have allocated
    //       enough space to do global operations (1 rank) but
    //       will fill a subset of that data that the rank owns.

    // Populate FABs from bdy interpolation (primitive vars)
    //==========================================================
    for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
        if (ivar == ivarQV && !l_relax_qv) { continue; }
        int ivar_idx = var_map[ivar];
        Box domain   = geom.Domain();
        auto ixtype  = S_cur_data[ivar_idx].boxArray().ixType();
        domain.convert(ixtype);
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // BndryReg idx and limiting (Rho has no bdy-plane slot)
        int bdy_comp = bnd_map[(ivar < RealBdyVars::NumTypes) ? ivar : int(RealBdyVars::U)];
        const auto& dom_cc_lo = lbound(geom.Domain());
        const auto& dom_cc_hi = ubound(geom.Domain());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_cur_data[ivar_idx],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // NOTE: Ghost cells needed for idx type mismatch between mask and data (do_upwind)
            IntVect ng_vect(0);
            //IntVect ng_vect(1,1,0);
            Box gtbx = grow(mfi.tilebox(ixtype.toIntVect()),ng_vect);
            Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
            realbdy_interior_bxs_xy(gtbx, domain, width,
                                    tbx_xlo, tbx_xhi,
                                    tbx_ylo, tbx_yhi,
                                    ng_vect, true);

            Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
            Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
            if (ivar  == ivarU) {
                arr_xlo = U_xlo.array(); arr_xhi = U_xhi.array();
                arr_ylo = U_ylo.array(); arr_yhi = U_yhi.array();
            } else if (ivar  == ivarV) {
                arr_xlo = V_xlo.array(); arr_xhi = V_xhi.array();
                arr_ylo = V_ylo.array(); arr_yhi = V_yhi.array();
            } else if (ivar  == ivarT){
                arr_xlo = T_xlo.array(); arr_xhi = T_xhi.array();
                arr_ylo = T_ylo.array(); arr_yhi = T_yhi.array();
            } else if (ivar  == ivarQV){
                arr_xlo = Q_xlo.array(); arr_xhi = Q_xhi.array();
                arr_ylo = Q_ylo.array(); arr_yhi = Q_yhi.array();
            } else if (ivar  == ivarR){
                arr_xlo = R_xlo.array(); arr_xhi = R_xhi.array();
                arr_ylo = R_ylo.array(); arr_yhi = R_yhi.array();
            } else {
                continue;
            }

            // Rho has no boundary plane -- its target is read from fcons_tgt
            // below. Bind the plane arrays to a valid slot so the (unused)
            // Array4s are well formed.
            const int pvar = (ivar < RealBdyVars::NumTypes) ? ivar : int(RealBdyVars::U);

            // Boundary data at fixed time intervals
            const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][pvar].const_array();
            const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][pvar].const_array();
            const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][pvar].const_array();
            const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][pvar].const_array();
            const auto& bdatylo_n   = bdy_data_ylo[n_time   ][pvar].const_array();
            const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][pvar].const_array();
            const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][pvar].const_array();
            const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][pvar].const_array();

            // Current density to convert to conserved vars
            Array4<Real> r_arr = S_cur_data[IntVars::cons].array(mfi);

            // Density used to couple the primitive bdy values, and the Rho
            // target itself. With the mass-consistent forcing both are the
            // ERA5 density rho*, so the momentum target is rho* u* and its
            // divergence is ERA5's rather than the model's. fcons_tgt shares
            // the cell-centered BoxArray/DM, so indexing it with this mfi is
            // the same pattern already used for r_arr above.
            Array4<const Real> rt_arr = (l_mass_consistent) ? fcons_tgt->const_array(mfi)
                                                            : Array4<const Real>{};
            const bool l_mc = l_mass_consistent;
            const int  l_iv = ivar;
            const int  l_ir = ivarR;
            // fcons_tgt's physical-boundary ghosts are zero-initialized, so
            // clamp into the valid region (the wall face then sees the same
            // zero-gradient density fill_from_realbdy imposes).
            const int cx_lo = dom_cc_lo.x, cx_hi = dom_cc_hi.x;
            const int cy_lo = dom_cc_lo.y, cy_hi = dom_cc_hi.y;

            // Limiting offset
            int offset = width - 1;

            // Populate with interpolation (protect from ghost cells)
            ParallelFor(tbx_xlo, tbx_xhi,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int ii = std::max(i , dom_lo.x); ii = std::min(ii, dom_lo.x+offset);
                int jj = std::max(j , dom_lo.y); jj = std::min(jj, dom_hi.y);

                // Coupling density: rho* under the mass-consistent forcing, so
                // the momentum target is rho* u* and its divergence is ERA5's
                // rather than carrying u* grad(rho_model) as a mass source.
                auto RHO = [=] (int a, int b, int c) -> Real {
                    if (!l_mc) { return r_arr(a,b,c); }
                    a = amrex::min(amrex::max(a,cx_lo),cx_hi);
                    b = amrex::min(amrex::max(b,cy_lo),cy_hi);
                    return rt_arr(a,b,c,Rho_comp);
                };
                Real rho_interp;
                if (l_iv==l_ir) {
                    rho_interp = RHO(i,j,k);            // Rho target is rho* itself
                } else if (ivar==ivarU) {
                    rho_interp = myhalf * ( RHO(i-1,j  ,k) + RHO(i,j,k) );
                } else if (ivar==ivarV) {
                    rho_interp = myhalf * ( RHO(i  ,j-1,k) + RHO(i,j,k) );
                } else {
                    rho_interp = RHO(i,j,k);
                }

                if (l_iv == l_ir) { arr_xlo(i,j,k) = rho_interp; return; }
                if (bdatxlo) {
                    int ii2 = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                    int jj2 = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                    arr_xlo(i,j,k) = rho_interp * bdatxlo(ii2,jj2,k,bdy_comp);
                } else {
                    arr_xlo(i,j,k) = rho_interp * ( oma   * bdatxlo_n  (ii,jj,k,0)
                                                  + alpha * bdatxlo_np1(ii,jj,k,0) );
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int ii = std::max(i , dom_hi.x-offset); ii = std::min(ii, dom_hi.x);
                int jj = std::max(j , dom_lo.y);        jj = std::min(jj, dom_hi.y);

                // Coupling density: rho* under the mass-consistent forcing, so
                // the momentum target is rho* u* and its divergence is ERA5's
                // rather than carrying u* grad(rho_model) as a mass source.
                auto RHO = [=] (int a, int b, int c) -> Real {
                    if (!l_mc) { return r_arr(a,b,c); }
                    a = amrex::min(amrex::max(a,cx_lo),cx_hi);
                    b = amrex::min(amrex::max(b,cy_lo),cy_hi);
                    return rt_arr(a,b,c,Rho_comp);
                };
                Real rho_interp;
                if (l_iv==l_ir) {
                    rho_interp = RHO(i,j,k);            // Rho target is rho* itself
                } else if (ivar==ivarU) {
                    rho_interp = myhalf * ( RHO(i-1,j  ,k) + RHO(i,j,k) );
                } else if (ivar==ivarV) {
                    rho_interp = myhalf * ( RHO(i  ,j-1,k) + RHO(i,j,k) );
                } else {
                    rho_interp = RHO(i,j,k);
                }

                if (l_iv == l_ir) { arr_xhi(i,j,k) = rho_interp; return; }
                if (bdatxhi) {
                    int ii2 = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                    int jj2 = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                    arr_xhi(i,j,k) = rho_interp * bdatxhi(ii2,jj2,k,bdy_comp);
                } else {
                    arr_xhi(i,j,k) = rho_interp * ( oma   * bdatxhi_n  (ii,jj,k,0)
                                                  + alpha * bdatxhi_np1(ii,jj,k,0) );
                }
            });

            ParallelFor(tbx_ylo, tbx_yhi,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int ii = std::max(i , dom_lo.x); ii = std::min(ii, dom_hi.x);
                int jj = std::max(j , dom_lo.y); jj = std::min(jj, dom_lo.y+offset);

                // Coupling density: rho* under the mass-consistent forcing, so
                // the momentum target is rho* u* and its divergence is ERA5's
                // rather than carrying u* grad(rho_model) as a mass source.
                auto RHO = [=] (int a, int b, int c) -> Real {
                    if (!l_mc) { return r_arr(a,b,c); }
                    a = amrex::min(amrex::max(a,cx_lo),cx_hi);
                    b = amrex::min(amrex::max(b,cy_lo),cy_hi);
                    return rt_arr(a,b,c,Rho_comp);
                };
                Real rho_interp;
                if (l_iv==l_ir) {
                    rho_interp = RHO(i,j,k);            // Rho target is rho* itself
                } else if (ivar==ivarU) {
                    rho_interp = myhalf * ( RHO(i-1,j  ,k) + RHO(i,j,k) );
                } else if (ivar==ivarV) {
                    rho_interp = myhalf * ( RHO(i  ,j-1,k) + RHO(i,j,k) );
                } else {
                    rho_interp = RHO(i,j,k);
                }

                if (l_iv == l_ir) { arr_ylo(i,j,k) = rho_interp; return; }
                if (bdatylo) {
                    int ii2 = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                    int jj2 = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                    arr_ylo(i,j,k) = rho_interp * bdatylo(ii2,jj2,k,bdy_comp);
                } else {
                    arr_ylo(i,j,k) = rho_interp * ( oma  * bdatylo_n  (ii,jj,k,0)
                                                  + alpha * bdatylo_np1(ii,jj,k,0) );
                }
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                int ii = std::max(i , dom_lo.x);        ii = std::min(ii, dom_hi.x);
                int jj = std::max(j , dom_hi.y-offset); jj = std::min(jj, dom_hi.y);

                // Coupling density: rho* under the mass-consistent forcing, so
                // the momentum target is rho* u* and its divergence is ERA5's
                // rather than carrying u* grad(rho_model) as a mass source.
                auto RHO = [=] (int a, int b, int c) -> Real {
                    if (!l_mc) { return r_arr(a,b,c); }
                    a = amrex::min(amrex::max(a,cx_lo),cx_hi);
                    b = amrex::min(amrex::max(b,cy_lo),cy_hi);
                    return rt_arr(a,b,c,Rho_comp);
                };
                Real rho_interp;
                if (l_iv==l_ir) {
                    rho_interp = RHO(i,j,k);            // Rho target is rho* itself
                } else if (ivar==ivarU) {
                    rho_interp = myhalf * ( RHO(i-1,j  ,k) + RHO(i,j,k) );
                } else if (ivar==ivarV) {
                    rho_interp = myhalf * ( RHO(i  ,j-1,k) + RHO(i,j,k) );
                } else {
                    rho_interp = RHO(i,j,k);
                }

                if (l_iv == l_ir) { arr_yhi(i,j,k) = rho_interp; return; }
                if (bdatyhi) {
                    int ii2 = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                    int jj2 = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                    arr_yhi(i,j,k) = rho_interp * bdatyhi(ii2,jj2,k,bdy_comp);
                } else {
                    arr_yhi(i,j,k) = rho_interp * ( oma   * bdatyhi_n  (ii,jj,k,0)
                                                  + alpha * bdatyhi_np1(ii,jj,k,0) );
                }
            });
        } // mfi
    } // ivar


    // Compute RHS in relaxation region
    //==========================================================
    auto dx = geom.CellSizeArray();
    auto ProbLo = geom.ProbLoArray();
    auto ProbHi = geom.ProbHiArray();

    for (int ivar(ivarU); ivar < BdyEnd; ivar++) {
        if (ivar == ivarQV && !l_relax_qv) { continue; }
        int ivar_idx = ivar_map[ivar];
        int icomp    = comp_map[ivar];

        Box domain = geom.Domain();
        domain.convert(S_cur_data[ivar_idx].boxArray().ixType());
        IntVect ng_vect(0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(S_cur_data[ivar_idx],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();
            Box tbx_xlo, tbx_xhi, tbx_ylo, tbx_yhi;
            realbdy_interior_bxs_xy(tbx, domain, width,
                                    tbx_xlo, tbx_xhi,
                                    tbx_ylo, tbx_yhi,
                                    ng_vect);

            Array4<Real> rhs_arr;  Array4<Real> data_arr;
            Array4<Real> arr_xlo;  Array4<Real> arr_xhi;
            Array4<Real> arr_ylo;  Array4<Real> arr_yhi;
            if (ivar  == ivarU) {
                arr_xlo  = U_xlo.array(); arr_xhi = U_xhi.array();
                arr_ylo  = U_ylo.array(); arr_yhi = U_yhi.array();
                rhs_arr  = S_rhs[IntVars::xmom].array(mfi);
                data_arr = S_cur_data[IntVars::xmom].array(mfi);
            } else if (ivar  == ivarV) {
                arr_xlo  = V_xlo.array(); arr_xhi = V_xhi.array();
                arr_ylo  = V_ylo.array(); arr_yhi = V_yhi.array();
                rhs_arr  = S_rhs[IntVars::ymom].array(mfi);
                data_arr = S_cur_data[IntVars::ymom].array(mfi);
            } else if (ivar  == ivarT){
                arr_xlo  = T_xlo.array(); arr_xhi = T_xhi.array();
                arr_ylo  = T_ylo.array(); arr_yhi = T_yhi.array();
                rhs_arr  = S_rhs[IntVars::cons].array(mfi);
                data_arr = S_cur_data[IntVars::cons].array(mfi);
            } else if (ivar  == ivarQV){
                arr_xlo  = Q_xlo.array(); arr_xhi = Q_xhi.array();
                arr_ylo  = Q_ylo.array(); arr_yhi = Q_yhi.array();
                rhs_arr  = S_rhs[IntVars::cons].array(mfi);
                data_arr = S_cur_data[IntVars::cons].array(mfi);
            } else if (ivar  == ivarR){
                // Rho relaxation toward rho*. This runs BEFORE the normal-face
                // product-rule overwrite below, so that overwrite picks up the
                // rho tendency including this term rather than a stale one.
                arr_xlo  = R_xlo.array(); arr_xhi = R_xhi.array();
                arr_ylo  = R_ylo.array(); arr_yhi = R_yhi.array();
                rhs_arr  = S_rhs[IntVars::cons].array(mfi);
                data_arr = S_cur_data[IntVars::cons].array(mfi);
            } else {
                continue;
            }

            realbdy_compute_relaxation(icomp, 1,
                                       width, dx, ProbLo, ProbHi, F1,
                                       tbx_xlo , tbx_xhi , tbx_ylo , tbx_yhi ,
                                       arr_xlo , arr_xhi , arr_ylo , arr_yhi ,
                                       data_arr, rhs_arr);
        } // mfi
    } // ivar

    // Cancel the relaxation ramp's spurious divergence
    //==========================================================
    // The relaxation adds S = F*F1*(A-B) to the MOMENTUM rhs, whose divergence
    // is F1*[ F div(A-B) + grad(F).(A-B) ]. The second term has no physical
    // counterpart: it is a divergence source proportional to the ramp gradient
    // times the local momentum error, and since rho is unconstrained the only
    // place it can go is Delta_z(rho Omega) -- i.e. it becomes w. Measured: it
    // peaks mid-band (d ~ 0.6*width) and its amplitude scales as 1/width.
    //
    // Cancel it with a companion MOMENTUM increment C obeying
    //     div C = -F1 * grad(F).(A-B)
    // so the rho equation is never given a source and remains a pure flux
    // divergence (this is what separates it from a volumetric rho relaxation).
    // grad(F) is wall-normal, so C is a 1-D wall-normal integral; with
    // F(d) = ((w-d)/w)^2 the discrete increment is dx-free:
    //     C(d+1) = C(d) + 2*F1*(w-d)/w^2 * err(d)
    // starting from C = 0 at the wall (the wall flux is left alone). The
    // integral has a nonzero net, which would otherwise leak momentum into the
    // interior, so a linear taper drives C back to zero at the inner edge. The
    // residual is then a UNIFORM divergence over the band whose net per column
    // is exactly what the barotropic wall flux correction absorbs -- the two
    // are complementary, not alternatives.
    //
    // Corners (within `width` of two walls) are left untreated: there F is
    // max(xi^2,eta^2) and the 1-D split does not apply.
    // erf.hindcast_ramp_div_correction: 0 = off,
    //   1 = TAPER   -- C=0 at both ends, residual smeared over the band.
    //                  MEASURED FAILURE: the net is weighted by (w-d)/w^2 times
    //                  an error that GROWS inward, so it peaks mid-band exactly
    //                  where the artifact does; the taper then adds back a
    //                  uniform divergence of comparable size. Peaks stayed
    //                  47.0/32.0/15.8 % at width 10/15/20 -- still ~1/width.
    //   2 = WALL    -- C=0 at the inner edge instead, so the residual leaves
    //                  through the wall face as a per-level mass flux. The
    //                  cancellation is then exact inside the band and nothing
    //                  leaks into the interior.
    static const int l_ramp_corr = [] {
        int m=0; amrex::ParmParse pp("erf");
        pp.query("hindcast_ramp_div_correction", m); return m; }();

    MultiFab cwx, cwy;
    if (l_ramp_corr == 2) {
        cwx.define(S_rhs[IntVars::xmom].boxArray(), S_rhs[IntVars::xmom].DistributionMap(), 1, 0);
        cwy.define(S_rhs[IntVars::ymom].boxArray(), S_rhs[IntVars::ymom].DistributionMap(), 1, 0);
        cwx.setVal(0.0); cwy.setVal(0.0);
    }

    if (l_ramp_corr > 0 && width > 1)
    {
        const auto rlo = lbound(geom.Domain());
        const auto rhi = ubound(geom.Domain());
        const int  wd  = width;
        const Real cF1 = F1;

        Array4<const Real> tgt_ulo = U_xlo.const_array();
        Array4<const Real> tgt_uhi = U_xhi.const_array();
        Array4<const Real> tgt_vlo = V_ylo.const_array();
        Array4<const Real> tgt_vhi = V_yhi.const_array();

        // ---- x walls ----
        for (MFIter mfi(S_rhs[IntVars::xmom]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            const Array4<Real>&       rhs = S_rhs[IntVars::xmom].array(mfi);
            const Array4<const Real>& mom = S_cur_data[IntVars::xmom].const_array(mfi);

            for (int side = 0; side < 2; ++side)
            {
                const bool lo = (side == 0);
                const int  iw = lo ? rlo.x : rhi.x + 1;          // wall face index
                if (lo  && vbx.smallEnd(0) != rlo.x)     { continue; }
                if (!lo && vbx.bigEnd(0)   != rhi.x + 1) { continue; }
                if (vbx.length(0) < wd + 1)              { continue; }

                Box slab = makeSlab(vbx, 0, iw);
                Array4<const Real> tgt = lo ? tgt_ulo : tgt_uhi;
                const int istep = lo ? 1 : -1;
                const int mode  = l_ramp_corr;
                Array4<Real> cw = (mode == 2) ? cwx.array(mfi) : Array4<Real>{};

                ParallelFor(slab, [=] AMREX_GPU_DEVICE (int, int j, int k) noexcept
                {
                    // Skip corners: the 1-D split is only valid where this wall
                    // is unambiguously the nearest one.
                    if (amrex::min(j - rlo.y, rhi.y - j) < wd) { return; }

                    Real net = Real(0.0);
                    for (int d = 0; d < wd; ++d) {
                        const int fa = iw + istep*d;
                        const int fb = iw + istep*(d+1);
                        const Real err = Real(0.5) * ( (tgt(fa,j,k) - mom(fa,j,k))
                                                     + (tgt(fb,j,k) - mom(fb,j,k)) );
                        net += Real(2.0) * cF1 * Real(wd-d) / Real(wd*wd) * err;
                    }

                    // mode 1: C(0)=0, taper to zero at the inner edge.
                    // mode 2: C(0)=-net so C(wd)=0 exactly; the residual exits
                    //         through the wall face instead of the band.
                    Real C = (mode == 2) ? -net : Real(0.0);
                    if (mode == 2) { cw(iw,j,k) = C; }
                    for (int d = 0; d <= wd; ++d) {
                        if (d > 0 || mode == 1) {
                            rhs(iw + istep*d, j, k) +=
                                (mode == 2) ? C : (C - net * Real(d) / Real(wd));
                        }
                        if (d < wd) {
                            const int fa = iw + istep*d;
                            const int fb = iw + istep*(d+1);
                            const Real err = Real(0.5) * ( (tgt(fa,j,k) - mom(fa,j,k))
                                                         + (tgt(fb,j,k) - mom(fb,j,k)) );
                            C += Real(2.0) * cF1 * Real(wd-d) / Real(wd*wd) * err;
                        }
                    }
                });
            }
        }

        // ---- y walls ----
        for (MFIter mfi(S_rhs[IntVars::ymom]); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.validbox();
            const Array4<Real>&       rhs = S_rhs[IntVars::ymom].array(mfi);
            const Array4<const Real>& mom = S_cur_data[IntVars::ymom].const_array(mfi);

            for (int side = 0; side < 2; ++side)
            {
                const bool lo = (side == 0);
                const int  jw = lo ? rlo.y : rhi.y + 1;
                if (lo  && vbx.smallEnd(1) != rlo.y)     { continue; }
                if (!lo && vbx.bigEnd(1)   != rhi.y + 1) { continue; }
                if (vbx.length(1) < wd + 1)              { continue; }

                Box slab = makeSlab(vbx, 1, jw);
                Array4<const Real> tgt = lo ? tgt_vlo : tgt_vhi;
                const int jstep = lo ? 1 : -1;
                const int mode  = l_ramp_corr;
                Array4<Real> cw = (mode == 2) ? cwy.array(mfi) : Array4<Real>{};

                ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int, int k) noexcept
                {
                    if (amrex::min(i - rlo.x, rhi.x - i) < wd) { return; }

                    Real net = Real(0.0);
                    for (int d = 0; d < wd; ++d) {
                        const int fa = jw + jstep*d;
                        const int fb = jw + jstep*(d+1);
                        const Real err = Real(0.5) * ( (tgt(i,fa,k) - mom(i,fa,k))
                                                     + (tgt(i,fb,k) - mom(i,fb,k)) );
                        net += Real(2.0) * cF1 * Real(wd-d) / Real(wd*wd) * err;
                    }

                    Real C = (mode == 2) ? -net : Real(0.0);
                    if (mode == 2) { cw(i,jw,k) = C; }
                    for (int d = 0; d <= wd; ++d) {
                        if (d > 0 || mode == 1) {
                            rhs(i, jw + jstep*d, k) +=
                                (mode == 2) ? C : (C - net * Real(d) / Real(wd));
                        }
                        if (d < wd) {
                            const int fa = jw + jstep*d;
                            const int fb = jw + jstep*(d+1);
                            const Real err = Real(0.5) * ( (tgt(i,fa,k) - mom(i,fa,k))
                                                         + (tgt(i,fb,k) - mom(i,fb,k)) );
                            C += Real(2.0) * cF1 * Real(wd-d) / Real(wd*wd) * err;
                        }
                    }
                });
            }
        }
    }

    // Set normal velocity RHS at the boundary
    //==========================================================
    // The barotropic wall flux correction is folded into u_val/v_val here as
    // well as into the Dirichlet value in fill_from_realbdy, so the product
    // rule keeps the wall momentum consistent with the density it is built
    // from instead of fighting the imposed value.
    static const bool l_wall_corr = [] {
        bool b=false; amrex::ParmParse pp("erf");
        pp.query("hindcast_wall_flux_correction", b); return b; }();
    static const Real l_wall_tau = [] {
        Real t=Real(3600.0); amrex::ParmParse pp("erf");
        pp.query("hindcast_wall_flux_tau", t); return t; }();
    const bool l_wc = l_wall_corr && (fcons_tgt != nullptr);

    MultiFab dvelc;
    if (l_wc) {
        dvelc.define(S_cur_data[IntVars::cons].boxArray(),
                     S_cur_data[IntVars::cons].DistributionMap(), 1, 1);
        compute_wall_flux_correction(geom, S_cur_data[IntVars::cons], *fcons_tgt,
                                     l_wall_tau, dvelc);
    }

    Box domain  = geom.Domain();
    Box domainx = convert(domain, IntVect(1,0,0));
    Box domainy = convert(domain, IntVect(0,1,0));

    int ilo = domainx.smallEnd(0);
    int ihi = domainx.bigEnd(0);
    int jlo = domainy.smallEnd(1);
    int jhi = domainy.bigEnd(1);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(S_cur_data[IntVars::cons],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box tbx = mfi.nodaltilebox(0);
        Box tbx_lo, tbx_hi;
        if (tbx.smallEnd(0) == ilo) {
            tbx_lo = makeSlab(tbx,0,ilo);
        }
        if (tbx.bigEnd(0) == ihi) {
            tbx_hi = makeSlab(tbx,0,ihi);
        }

        Box tby = mfi.nodaltilebox(1);
        Box tby_lo, tby_hi;
        if (tby.smallEnd(1) == jlo) {
            tby_lo = makeSlab(tby,1,jlo);
        }
        if (tby.bigEnd(1) == jhi) {
            tby_hi = makeSlab(tby,1,jhi);
        }

        Array4<Real> rhs_xmom  = S_rhs[IntVars::xmom].array(mfi);
        Array4<Real> rhs_ymom  = S_rhs[IntVars::ymom].array(mfi);

        Array4<const Real> rhs_cons = S_rhs[IntVars::cons].const_array(mfi);
        Array4<const Real> cons_arr = S_cur_data[IntVars::cons].const_array(mfi);
        Array4<const Real> dvarr    = l_wc ? dvelc.const_array(mfi) : Array4<const Real>{};
        const auto wlo = lbound(geom.Domain());
        const auto whi = ubound(geom.Domain());
        const bool l_wc_l = l_wc;
        // Ramp-cancellation residual routed out through the wall (mode 2).
        // This block ASSIGNS rhs at the wall faces, so the residual has to be
        // added here rather than left for the correction pass above.
        const bool l_rc2 = (l_ramp_corr == 2);
        Array4<const Real> cwx_a = l_rc2 ? cwx.const_array(mfi) : Array4<const Real>{};
        Array4<const Real> cwy_a = l_rc2 ? cwy.const_array(mfi) : Array4<const Real>{};

        const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][ivarU].const_array();
        const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][ivarU].const_array();
        const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][ivarU].const_array();
        const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][ivarU].const_array();

        const auto& bdatylo_n   = bdy_data_ylo[n_time   ][ivarV].const_array();
        const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][ivarV].const_array();
        const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][ivarV].const_array();
        const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][ivarV].const_array();

        ParallelFor(tbx_lo, tbx_hi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho_tend = rhs_cons(i,j,k);
            Real rho_val  = Real(0.5) * (cons_arr(i,j,k) + cons_arr(i-1,j,k));
            Real u_tend, u_val;
            if (btenxlo) {
                u_tend   = btenxlo(i,j,k,BCVars::xvel_bc);
                u_val    = bdatxlo(i,j,k,BCVars::xvel_bc);
            } else {
                u_tend   = (bdatxlo_np1(i,j,k) - bdatxlo_n(i,j,k)) / bdy_time_interval;
                u_val    = oma * bdatxlo_n(i,j,k) + alpha * bdatxlo_np1(i,j,k);
            }
            if (l_wc_l) { u_val += dvarr(wlo.x, amrex::min(amrex::max(j,wlo.y),whi.y), k); }
            rhs_xmom(i,j,k) = rho_val * u_tend + u_val * rho_tend;
            if (l_rc2) { rhs_xmom(i,j,k) += cwx_a(i,j,k); }
            },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // i == ihi is the x-face at nx, so cell (i,j,k) is the FIRST GHOST:
            // S_rhs[cons] there is never written (every RHS writer works on
            // mfi.tilebox()) and F_slow is setVal(0) once at construction, so
            // reading rhs_cons(i,...) silently drops the u*drho/dt term on the
            // hi walls while it stays live on the lo walls. Use the adjacent
            // valid cell, mirroring the lo side.
            Real rho_tend = rhs_cons(i-1,j,k);
            Real rho_val  = Real(0.5) * (cons_arr(i,j,k) + cons_arr(i-1,j,k));
            Real u_tend, u_val;
            if (btenxhi) {
                u_tend   = btenxhi(i,j,k,BCVars::xvel_bc);
                u_val    = bdatxhi(i,j,k,BCVars::xvel_bc);
            } else {
                u_tend   = (bdatxhi_np1(i,j,k) - bdatxhi_n(i,j,k)) / bdy_time_interval;
                u_val    = oma * bdatxhi_n(i,j,k) + alpha * bdatxhi_np1(i,j,k);
            }
            if (l_wc_l) { u_val -= dvarr(whi.x, amrex::min(amrex::max(j,wlo.y),whi.y), k); }
            rhs_xmom(i,j,k) = rho_val * u_tend + u_val * rho_tend;
            if (l_rc2) { rhs_xmom(i,j,k) += cwx_a(i,j,k); }
        });

        ParallelFor(tby_lo, tby_hi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho_tend = rhs_cons(i,j,k);
            Real rho_val  = Real(0.5) * (cons_arr(i,j,k) + cons_arr(i,j-1,k));
            Real v_tend, v_val;
            if (btenylo) {
                v_tend   = btenylo(i,j,k,BCVars::yvel_bc);
                v_val    = bdatylo(i,j,k,BCVars::yvel_bc);
            } else {
                v_tend   = (bdatylo_np1(i,j,k) - bdatylo_n(i,j,k)) / bdy_time_interval;
                v_val    = oma * bdatylo_n(i,j,k) + alpha * bdatylo_np1(i,j,k);
            }
            if (l_wc_l) { v_val += dvarr(amrex::min(amrex::max(i,wlo.x),whi.x), wlo.y, k); }
            rhs_ymom(i,j,k) = rho_val * v_tend + v_val * rho_tend;
            if (l_rc2) { rhs_ymom(i,j,k) += cwy_a(i,j,k); }
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // See the x-hi note: (i,j,k) here is the first ghost in y.
            Real rho_tend = rhs_cons(i,j-1,k);
            Real rho_val  = Real(0.5) * (cons_arr(i,j,k) + cons_arr(i,j-1,k));
            Real v_tend, v_val;
            if (btenyhi) {
                v_tend   = btenyhi(i,j,k,BCVars::yvel_bc);
                v_val    = bdatyhi(i,j,k,BCVars::yvel_bc);
            } else {
                v_tend   = (bdatyhi_np1(i,j,k) - bdatyhi_n(i,j,k)) / bdy_time_interval;
                v_val    = oma * bdatyhi_n(i,j,k) + alpha * bdatyhi_np1(i,j,k);
            }
            if (l_wc_l) { v_val -= dvarr(amrex::min(amrex::max(i,wlo.x),whi.x), whi.y, k); }
            rhs_ymom(i,j,k) = rho_val * v_tend + v_val * rho_tend;
            if (l_rc2) { rhs_ymom(i,j,k) += cwy_a(i,j,k); }
        });
    } // mfi
}

/**
 * Compute the RHS in the fine relaxation zone
 *
 * @param[in]  time      current (elapsed) time
 * @param[in]  delta_t   timestep
 * @param[in]  width     number of cells in (relaxation+specified) zone
 * @param[in]  set_width number of cells in (specified) zone
 * @param[in]  FPr_c     cons fine patch container
 * @param[in]  FPr_u     uvel fine patch container
 * @param[in]  FPr_v     vvel fine patch container
 * @param[in]  FPr_w     wvel fine patch container
 * @param[in]  boxes_at_level boxes at current level
 * @param[in]  domain_bcs_type boundary condition types
 * @param[out] S_rhs     RHS to be computed here
 * @param[in]  S_data    current value of the solution
 */
void
fine_compute_interior_ghost_rhs (const Real& time,
                                 const Real& delta_t,
                                 const int& width,
                                 const int& set_width,
                                 const Geometry& geom,
                                 ERFFillPatcher* FPr_c,
                                 ERFFillPatcher* FPr_u,
                                 ERFFillPatcher* FPr_v,
                                 ERFFillPatcher* FPr_w,
                                 Vector<BCRec>& domain_bcs_type,
                                 Vector<MultiFab>& S_rhs_f,
                                 Vector<MultiFab>& S_data_f)
{
    BL_PROFILE_REGION("fine_compute_interior_ghost_RHS()");

    // Relaxation constants. The upstream (dead-code) constants gave a
    // tau = 10*dt momentum relaxation -- strong enough that under a winter
    // jet the band's forced momenta violate discrete continuity faster than
    // the density blend can absorb (measured: band rho 1.2 -> 2.5 in ~2.5
    // min, Jan-9). erf.cf_nudge_factor relaxes that timescale.
    static const Real nudge_fac = [] {
        Real f = Real(10.); amrex::ParmParse pp("erf"); pp.query("cf_nudge_factor", f); return f;
    }();
    Real F1 = one/(nudge_fac*delta_t);
    Real F2 = one/(Real(5.)*nudge_fac*delta_t);

    // Vector of MFs to hold data (dm differs w/ fine patch)
    Vector<MultiFab> fmf_p_v;

    // Loop over the variables
    for (int ivar_idx = 0; ivar_idx < IntVars::NumTypes; ++ivar_idx)
    {
        // Fine mfs
        MultiFab& fmf = S_data_f[ivar_idx];
        MultiFab& rhs = S_rhs_f [ivar_idx];

        // NOTE: These temporary MFs and copy operations are horrible
        //       for memory usage and efficiency. However, we need to
        //       have access to ghost cells in the cons array to convert
        //       from primitive u/v/w to momentum. Furthermore, the BA
        //       for the fine patches in ERFFillPatcher don't match the
        //       BA for the data/RHS. For this reason, the data is copied
        //       to a vector of MFs (with ghost cells) so the BAs match
        //       the BA of data/RHS and we have access to rho to convert
        //       prim to conserved.

        // Temp MF on box (distribution map differs w/ fine patch)
        int num_var = fmf.nComp();
        fmf_p_v.emplace_back(fmf.boxArray(), fmf.DistributionMap(), num_var, fmf.nGrowVect());
        MultiFab& fmf_p = fmf_p_v[ivar_idx];
        MultiFab::Copy(fmf_p,fmf, 0, 0, num_var, fmf.nGrowVect());

        // Integer mask MF
        int set_mask_val;
        int relax_mask_val;
        iMultiFab* mask;

        // Fill fine patch on interior halo region
        //==========================================================
        if (ivar_idx == IntVars::cons)
        {
            FPr_c->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_c->GetMask();
            set_mask_val   = FPr_c->GetSetMaskVal();
            relax_mask_val = FPr_c->GetRelaxMaskVal();
        }
        // NOTE: RegisterCoarseData in ERF_Advance.cpp registers the coarse
        // MOMENTA (state_old/new[IntVars::xmom] etc.), so the FillRelax
        // result here is already momentum. The original (dead) version of
        // this routine multiplied by an interpolated rho, which would give
        // rho^2 * u -- that conversion is removed.
        else if (ivar_idx == IntVars::xmom)
        {
            FPr_u->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_u->GetMask();
            set_mask_val   = FPr_u->GetSetMaskVal();
            relax_mask_val = FPr_u->GetRelaxMaskVal();
        }
        else if (ivar_idx == IntVars::ymom)
        {
            FPr_v->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_v->GetMask();
            set_mask_val   = FPr_v->GetSetMaskVal();
            relax_mask_val = FPr_v->GetRelaxMaskVal();
        }
        else if (ivar_idx == IntVars::zmom)
        {
            FPr_w->FillRelax(fmf_p, time, void_bc, domain_bcs_type);
            mask           = FPr_w->GetMask();
            set_mask_val   = FPr_w->GetSetMaskVal();
            relax_mask_val = FPr_w->GetRelaxMaskVal();
        } else {
            Abort("Dont recognize this variable type in fine_compute_interior_ghost_RHS");
        }


        // Zero RHS in set region
        //==========================================================
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(rhs,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();
            const Array4<Real>& rhs_arr  = rhs.array(mfi);
            const Array4<const int>& mask_arr = mask->const_array(mfi);

            ParallelFor(tbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (mask_arr(i,j,k) == set_mask_val) {
                    rhs_arr(i,j,k) = zero;
                }
            });
        } // mfi

        // For Laplacian stencil
        rhs.FillBoundary(geom.periodicity());


        // Compute RHS in relaxation region
        //==========================================================
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(fmf_p,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box tbx = mfi.tilebox();
            const Array4<Real>&        rhs_arr = rhs.array(mfi);
            const Array4<const Real>& fine_arr = fmf_p.const_array(mfi);
            const Array4<const Real>& data_arr = fmf.const_array(mfi);
            const Array4<const int>&  mask_arr = mask->const_array(mfi);

            Box vbx = mfi.validbox();
            const auto& vbx_lo = lbound(vbx);
            const auto& vbx_hi = ubound(vbx);

            int icomp = 0;

            int Spec_z  = set_width;
            int Relax_z = width - Spec_z;
            Real num    = Real(Spec_z + Relax_z);
            Real denom  = Real(Relax_z - 1);
            // erf.cf_relax_rho=false skips relaxing Rho through the slow RHS
            // (a mass-equation relaxation can destabilize the acoustic
            // substepping under strong cross-boundary flow -- measured:
            // nest SE-corner rho blowup ~150 fine steps, Jan-9 jet).
            // DEFAULT true: rho relaxation is part of the configuration
            // verified stable on the Sep-9 convective case (3,500 steps
            // zero-warning); skipping it lets band mass ratchet instead.
            static const bool relax_rho = [] {
                bool b = true; amrex::ParmParse pp("erf"); pp.query("cf_relax_rho", b); return b;
            }();
            const bool skip_rho = (!relax_rho) && (ivar_idx == IntVars::cons);
            ParallelFor(tbx, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               if (skip_rho && (n + icomp == Rho_comp)) { return; }
               if (mask_arr(i,j,k) == relax_mask_val) {

                   // Indices
                   Real n_ind(-1); // Set to -1 to quiet compiler warning
                   int ii(width-1); int jj(width-1);
                   bool near_x_lo_wall(false); bool near_x_hi_wall(false);
                   bool near_y_lo_wall(false); bool near_y_hi_wall(false);
                   bool mask_x_found(false);   bool mask_y_found(false);

                   // Near x-wall
                   if ((i-vbx_lo.x) < width) {
                       near_x_lo_wall = true;
                       ii = i-vbx_lo.x;
                       if (mask_arr(vbx_lo.x,j,k) == 2) mask_x_found = true;
                   } else if ((vbx_hi.x-i) < width) {
                       near_x_hi_wall = true;
                       ii = vbx_hi.x-i;
                       if (mask_arr(vbx_hi.x,j,k) == 2) mask_x_found = true;
                   }

                   // Near y-wall
                   if ((j-vbx_lo.y) < width) {
                       near_y_lo_wall = true;
                       jj = j-vbx_lo.y;
                       if (mask_arr(i,vbx_lo.y,k) == 2) mask_y_found = true;
                   } else if ((vbx_hi.y-j) < width) {
                       near_y_hi_wall = true;
                       jj = vbx_hi.y-j;
                       if (mask_arr(i,vbx_hi.y,k) == 2) mask_y_found = true;
                   }

                   // Found a nearby masked cell (valid n_ind)
                   if (mask_x_found && mask_y_found) {
                       n_ind = std::min(ii,jj) + one;
                   } else if (mask_x_found) {
                       n_ind = ii + one;
                   } else if (mask_y_found) {
                       n_ind = jj + one;
                   // Pesky corner cell
                   } else {
                       if (near_x_lo_wall || near_x_hi_wall) {
                           Real dj_min{width-one};
                           int j_lb = std::max(vbx_lo.y,j-width);
                           int j_ub = std::min(vbx_hi.y,j+width);
                           int li   = (near_x_lo_wall) ? vbx_lo.x : vbx_hi.x;
                           for (int lj(j_lb); lj<=j_ub; ++lj) {
                               if (mask_arr(li,lj,k) == 2) {
                                   mask_y_found = true;
                                   dj_min = std::min(dj_min,(Real) std::abs(lj-j));
                               }
                           }
                           if (mask_y_found) {
                               Real mag = std::sqrt( Real(dj_min*dj_min + ii*ii) );
                               n_ind = std::min(mag,width-one) + one;
                           } else {
                               Abort("Mask not found near x wall!");
                           }
                       } else if (near_y_lo_wall || near_y_hi_wall) {
                           Real di_min{width-one};
                           int i_lb = std::max(vbx_lo.x,i-width);
                           int i_ub = std::min(vbx_hi.x,i+width);
                           int lj   = (near_y_lo_wall) ? vbx_lo.y : vbx_hi.y;
                           for (int li(i_lb); li<=i_ub; ++li) {
                               if (mask_arr(li,lj,k) == 2) {
                                   mask_x_found = true;
                                   di_min = std::min(di_min,(Real) std::abs(li-i));
                               }
                           }
                           if (mask_x_found) {
                               Real mag = std::sqrt( Real(di_min*di_min + jj*jj) );
                               n_ind = std::min(mag,width-one) + one;
                           } else {
                               Abort("Mask not found near y wall!");
                           }
                       } else {
                           Abort("Relaxation cell must be near a wall!");
                       }
                   }

                   Real Factor   = (num - n_ind)/denom;
                   Real d        = data_arr(i  ,j  ,k  ,n+icomp) + delta_t*rhs_arr(i  , j  , k  ,n+icomp);
                   Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i+1, j  , k  ,n+icomp);
                   Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i-1, j  , k  ,n+icomp);
                   Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp) + delta_t*rhs_arr(i  , j+1, k  ,n+icomp);
                   Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp) + delta_t*rhs_arr(i  , j-1, k  ,n+icomp);
                   Real delta    = fine_arr(i  ,j  ,k,n) - d;
                   Real delta_xp = fine_arr(i+1,j  ,k,n) - d_ip1;
                   Real delta_xm = fine_arr(i-1,j  ,k,n) - d_im1;
                   Real delta_yp = fine_arr(i  ,j+1,k,n) - d_jp1;
                   Real delta_ym = fine_arr(i  ,j-1,k,n) - d_jm1;
                   Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - Real(4.0)*delta;
                   rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;
               }
            });
        } // mfi
    } // ivar_idx
}
