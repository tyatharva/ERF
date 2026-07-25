#include "ERF.H"
#include "ERF_Utils.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF
/*
 * Impose boundary conditions using data read in from wrfbdy files
 *
 * @param[out] mfs  Vector of MultiFabs to be filled
 * @param[in] time  time at which the data should be filled
 */

void
ERF::fill_from_realbdy (const Vector<MultiFab*>& mfs,
                        const Real time,
                        bool cons_only,
                        int icomp_cons,
                        int ncomp_cons,
                        IntVect ngvect_cons,
                        IntVect ngvect_vels)
{
    int lev = 0;

    // We do not operate on the z ghost cells
    ngvect_cons[2] = 0;
    ngvect_vels[2] = 0;

    // Time interpolation
    Real dT = bdy_time_interval;

    Real time_tot = time + start_time;
    Real time_since_start_bdy = time_tot - start_bdy_time;
    int n_time    = static_cast<int>( time_since_start_bdy /  dT);
    int n_time_p1 = n_time + 1;
    Real alpha    = (time_since_start_bdy - n_time * dT) / dT;

    // Do not over run the last bdy file
    if (time_tot >= final_bdy_time) {
        n_time    = static_cast<int>( (final_bdy_time - start_bdy_time)/ dT);
        n_time_p1 = n_time;
        alpha     = zero;
    }

    AMREX_ALWAYS_ASSERT( alpha >= zero && alpha <= one);
    Real oma   = one - alpha;

    // Flags for read vars and index mapping
    Vector<int> cons_read = {0, 1, 0, 0,
                             1, 0, 0,
                             0, 0, 0,
                             0, 0, 0,
                             0, 0};
    // HINDCAST ONLY (UPSTREAM_ISSUES #14): also specify rho at the wall from
    // the hindcast boundary planes, so the imposed wall mass flux is rho* u*
    // -- the same flux the band interior is relaxed toward -- rather than
    // rho_model u* built on a zero-gradient density. Gated on a knob no
    // metgrid/wrfbdy deck sets, because HindcastBdyVars::RHO deliberately
    // aliases WRFBdyVars::PH. Implied by the mass-consistent lateral forcing,
    // still separately settable for A/B.
    //
    // w is deliberately NOT specifiable here: the frame's w slot carries ERA5
    // omega in Pa/s, not a geometric vertical velocity, and the correct Omega
    // follows from the mass budget closing rather than from imposition (cf.
    // the abandoned lateral WfromOmega experiment, upstream PR #2872).
    static const bool l_bdy_rho = [] {
        bool b=false, mc=false; ParmParse pp("erf");
        pp.query("hindcast_bdy_rho", b);
        pp.query("hindcast_mass_consistent_bdy", mc);
        return b || mc; }();

    // MPAS-A precedent: set w = 0 in the specified zone instead of taking a
    // zero-gradient copy from the interior. Reported to alleviate spurious
    // streamers and near-boundary instability under strong inflow.
    static const bool l_w_zero = [] {
        bool b=false; ParmParse pp("erf");
        pp.query("hindcast_bdy_w_zero", b); return b; }();
    const bool l_have_ext = (!bdy_data_xlo.empty()) &&
        (static_cast<int>(bdy_data_xlo[0].size()) > HindcastBdyVars::RHO);
    if (l_bdy_rho && l_have_ext) { cons_read[Rho_comp] = 1; }

    Vector<Vector<int>> is_read;
    is_read.push_back( cons_read );
    is_read.push_back( {1} ); // xvel
    is_read.push_back( {1} ); // yvel
    is_read.push_back( {0} ); // zvel

    // Real BC mapping (WRF/MetGrid)
    Vector<int> cons_map = {HindcastBdyVars::RHO, RealBdyVars::T, RhoKE_comp, RhoScalar_comp,
                            RealBdyVars::QV, RhoQ2_comp, RhoQ3_comp,
                            RhoQ4_comp, RhoQ5_comp, RhoQ6_comp,
                            RhoQ7_comp, RhoQ8_comp, RhoQ9_comp,
                            RhoQ10_comp, RhoQ11_comp};
    Vector<Vector<int>> ind_map;
    ind_map.push_back( cons_map );
    ind_map.push_back( {RealBdyVars::U} ); // xvel
    ind_map.push_back( {RealBdyVars::V} ); // yvel
    ind_map.push_back( {0} );                  // zvel (never read; see note above)

    // Bndry plane mapping
    Vector<int> bnd_cons_map = {Rho_comp, BCVars::RhoTheta_bc_comp, RhoKE_comp, RhoScalar_comp,
                                BCVars::RhoQ1_bc_comp, RhoQ2_comp, RhoQ3_comp,
                                RhoQ4_comp, RhoQ5_comp, RhoQ6_comp,
                                RhoQ7_comp, RhoQ8_comp, RhoQ9_comp,
                                RhoQ10_comp, RhoQ11_comp};
    Vector<Vector<int>> bnd_ind_map;
    bnd_ind_map.push_back( bnd_cons_map );
    bnd_ind_map.push_back( {BCVars::xvel_bc} ); // xvel
    bnd_ind_map.push_back( {BCVars::yvel_bc} ); // yvel
    bnd_ind_map.push_back( {0} );               // zvel
    Array4<Real> bdatxlo, bdatxhi, bdatylo, bdatyhi;
    if (m_r2d) {
        // Index is [plane orientation] and [level]
        Vector<std::unique_ptr<PlaneVector>>& bndry_data = m_r2d->interp_in_time(time_tot);
        bdatxlo = (*bndry_data[0])[0].array();
        bdatylo = (*bndry_data[1])[0].array();
        bdatxhi = (*bndry_data[3])[0].array();
        bdatyhi = (*bndry_data[4])[0].array();
    }

    // Nvars to loop over
    Vector<int> comp_var = {ncomp_cons, 1, 1, 1};

    // End of vars loop
    int var_idx_end = (cons_only) ? Vars::cons + 1 : Vars::NumTypes;

    // Loop over all variable types
    for (int var_idx = Vars::cons; var_idx < var_idx_end; ++var_idx)
    {
        MultiFab& mf = *mfs[var_idx];

        mf.FillBoundary(geom[lev].periodicity());

        // Note that "domain" is mapped onto the type of box the data is in
        Box domain = geom[lev].Domain();
        domain.convert(mf.boxArray().ixType());
        const auto& dom_lo = lbound(domain);
        const auto& dom_hi = ubound(domain);

        // BndryReg idx limiting
        const auto& dom_cc_lo = lbound(geom[lev].Domain());
        const auto& dom_cc_hi = ubound(geom[lev].Domain());

        // Offset only applies to cons (we may fill a subset of these vars)
        int offset = (var_idx == Vars::cons) ? icomp_cons : 0;

        // Ghost cells to be filled
        IntVect ng_vect = (var_idx == Vars::cons) ? ngvect_cons : ngvect_vels;

        // Set region width
        int set_width = 1;

        // Loop over each component
        for (int comp_idx(offset); comp_idx < (comp_var[var_idx]+offset); ++comp_idx)
        {

            // Variable can be read from wrf bdy
            //------------------------------------
            if (is_read[var_idx][comp_idx])
            {
                int ivar    = ind_map[var_idx][comp_idx];
                int bnd_var = bnd_ind_map[var_idx][comp_idx];

                // We have data at fixed time intervals we will call dT
                // Then to interpolate, given time, we can define n = (time/dT)
                // and alpha = (time - n*dT) / dT, then we define the data at time
                // as  alpha * (data at time n+1) + (1 - alpha) * (data at time n)
                const auto& bdatxlo_n   = bdy_data_xlo[n_time   ][ivar].const_array();
                const auto& bdatxlo_np1 = bdy_data_xlo[n_time_p1][ivar].const_array();
                const auto& bdatxhi_n   = bdy_data_xhi[n_time   ][ivar].const_array();
                const auto& bdatxhi_np1 = bdy_data_xhi[n_time_p1][ivar].const_array();
                const auto& bdatylo_n   = bdy_data_ylo[n_time   ][ivar].const_array();
                const auto& bdatylo_np1 = bdy_data_ylo[n_time_p1][ivar].const_array();
                const auto& bdatyhi_n   = bdy_data_yhi[n_time   ][ivar].const_array();
                const auto& bdatyhi_np1 = bdy_data_yhi[n_time_p1][ivar].const_array();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Grown tilebox so we fill exterior ghost cells as well
                    Box gbx = mfi.growntilebox(ng_vect);
                    const Array4<Real>& dest_arr = mf.array(mfi);
                    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                    realbdy_bc_bxs_xy(gbx, domain, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect);

                    // x-faces (includes exterior y ghost cells)
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (bdatxlo) {
                            int ii = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                            int jj = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                            dest_arr(i,j,k,comp_idx) = bdatxlo(ii,jj,k,bnd_var);
                        } else {
                            int ii = std::max(i , dom_lo.x);
                            int jj = std::max(j , dom_lo.y);
                                jj = std::min(jj, dom_hi.y);
                                dest_arr(i,j,k,comp_idx) = oma   * bdatxlo_n  (ii,jj,k,0)
                                                         + alpha * bdatxlo_np1(ii,jj,k,0);
                        }
                        // Rho itself is stored directly; everything else in cons
                        // is a rho-weighted quantity built from a plain bdy value.
                        if (var_idx == Vars::cons && comp_idx != Rho_comp) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (bdatxhi) {
                            int ii = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                            int jj = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                            dest_arr(i,j,k,comp_idx) = bdatxhi(ii,jj,k,bnd_var);
                        } else {
                            int ii = std::min(i , dom_hi.x);
                            int jj = std::max(j , dom_lo.y);
                                jj = std::min(jj, dom_hi.y);
                                dest_arr(i,j,k,comp_idx) = oma   * bdatxhi_n  (ii,jj,k,0)
                                                         + alpha * bdatxhi_np1(ii,jj,k,0);
                        }
                        // Rho itself is stored directly; everything else in cons
                        // is a rho-weighted quantity built from a plain bdy value.
                        if (var_idx == Vars::cons && comp_idx != Rho_comp) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    });

                    // y-faces (do not include exterior x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (bdatylo) {
                            int ii = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                            int jj = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                            dest_arr(i,j,k,comp_idx) = bdatylo(ii,jj,k,bnd_var);
                        } else {
                            int jj = std::max(j , dom_lo.y);
                            dest_arr(i,j,k,comp_idx) = oma   * bdatylo_n  (i,jj,k,0)
                                                     + alpha * bdatylo_np1(i,jj,k,0);
                        }
                        // Rho itself is stored directly; everything else in cons
                        // is a rho-weighted quantity built from a plain bdy value.
                        if (var_idx == Vars::cons && comp_idx != Rho_comp) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        if (bdatyhi) {
                            int ii = std::min(std::max(i , dom_cc_lo.x), dom_cc_hi.x);
                            int jj = std::min(std::max(j , dom_cc_lo.y), dom_cc_hi.y);
                            dest_arr(i,j,k,comp_idx) = bdatyhi(ii,jj,k,bnd_var);
                        } else {
                            int jj = std::min(j , dom_hi.y);
                            dest_arr(i,j,k,comp_idx) = oma   * bdatyhi_n  (i,jj,k,0)
                                                     + alpha * bdatyhi_np1(i,jj,k,0);
                        }
                        // Rho itself is stored directly; everything else in cons
                        // is a rho-weighted quantity built from a plain bdy value.
                        if (var_idx == Vars::cons && comp_idx != Rho_comp) dest_arr(i,j,k,comp_idx) *= dest_arr(i,j,k,Rho_comp);
                    });
                } // mfi

            // Variable not read from wrf bdy
            //------------------------------------
            } else {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    // Grown tilebox so we fill exterior ghost cells as well
                    Box gbx = mfi.growntilebox(ng_vect);
                    Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
                    realbdy_bc_bxs_xy(gbx, domain, set_width,
                                      bx_xlo, bx_xhi,
                                      bx_ylo, bx_yhi,
                                      ng_vect);

                    // Bounding
                    int i_xlo = bx_xlo.bigEnd(0)   + set_width;
                    int i_xhi = bx_xhi.smallEnd(0) - set_width;
                    int j_ylo = bx_ylo.bigEnd(1)   + set_width;
                    int j_yhi = bx_yhi.smallEnd(1) - set_width;

                    // Destination array
                    const Array4<Real>& dest_arr = mf.array(mfi);

                    // Moisture beyond qv is zeroed in the specified zone; with
                    // erf.hindcast_bdy_w_zero, w is too (MPAS-A precedent).
                    const bool zero_here = (comp_idx > RhoQ1_comp) ||
                                           (l_w_zero && (var_idx == Vars::zvel));

                    // x-faces (includes y ghost cells)
                    ParallelFor(bx_xlo, bx_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = zero_here ? Real(0.) :
                                                                            dest_arr(i_xlo,jj,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        int jj = std::max(j , dom_lo.y);
                            jj = std::min(jj, dom_hi.y);
                        dest_arr(i,j,k,comp_idx) = zero_here ? Real(0.) :
                                                                             dest_arr(i_xhi,jj,k,comp_idx);
                    });

                    // y-faces (does not include x ghost cells)
                    ParallelFor(bx_ylo, bx_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = zero_here ? Real(0.) :
                                                                             dest_arr(i,j_ylo,k,comp_idx);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        dest_arr(i,j,k,comp_idx) = zero_here ? Real(0.) :
                                                                             dest_arr(i,j_yhi,k,comp_idx);
                    });
                } // mfi
            } // is_read
        } // comp
    } // var

    // Barotropic wall mass-flux correction (erf.hindcast_wall_flux_correction).
    //
    // Delivers the column-mass constraint through the wall FLUX rather than a
    // volumetric source in the rho equation: the wall-normal velocity picks up
    // a single depth-independent increment sized so the imposed flux carries
    // the column's mass error. The rho equation keeps a pure flux divergence,
    // so the acoustic solver is never forced, and the increment has no
    // vertical structure to project onto vertical modes.
    static const bool l_wall_corr = [] {
        bool b=false; ParmParse pp("erf");
        pp.query("hindcast_wall_flux_correction", b); return b; }();
    static const Real l_wall_tau = [] {
        Real t=Real(3600.0); ParmParse pp("erf");
        pp.query("hindcast_wall_flux_tau", t); return t; }();

    if (l_wall_corr && !cons_only &&
        solverChoice.init_type == InitType::HindCast &&
        !forecast_state_interp[lev].empty())
    {
        const MultiFab& cons_mf = *mfs[Vars::cons];
        MultiFab dvel(cons_mf.boxArray(), cons_mf.DistributionMap(), 1, 1);
        compute_wall_flux_correction(geom[lev], cons_mf,
                                     forecast_state_interp[lev][Vars::cons],
                                     l_wall_tau, dvel);

        const auto dlo = lbound(geom[lev].Domain());
        const auto dhi = ubound(geom[lev].Domain());
        const IntVect ngv = ngvect_vels;

        for (int vdir = 0; vdir < 2; ++vdir)
        {
            MultiFab& mf = (vdir == 0) ? *mfs[Vars::xvel] : *mfs[Vars::yvel];
            Box domv = geom[lev].Domain();
            domv.convert(mf.boxArray().ixType());

            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(ngv);
                Box b_xlo, b_xhi, b_ylo, b_yhi;
                realbdy_bc_bxs_xy(gbx, domv, 1, b_xlo, b_xhi, b_ylo, b_yhi, ngv);

                const Array4<Real>&       vel = mf.array(mfi);
                const Array4<const Real>& dv  = dvel.const_array(mfi);

                if (vdir == 0) {
                    // u is the normal velocity at the x walls; + is inward at
                    // x-lo and outward at x-hi.
                    ParallelFor(b_xlo, b_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        int jj = amrex::min(amrex::max(j,dlo.y),dhi.y);
                        vel(i,j,k) += dv(dlo.x,jj,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        int jj = amrex::min(amrex::max(j,dlo.y),dhi.y);
                        vel(i,j,k) -= dv(dhi.x,jj,k);
                    });
                } else {
                    ParallelFor(b_ylo, b_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        int ii = amrex::min(amrex::max(i,dlo.x),dhi.x);
                        vel(i,j,k) += dv(ii,dlo.y,k);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        int ii = amrex::min(amrex::max(i,dlo.x),dhi.x);
                        vel(i,j,k) -= dv(ii,dhi.y,k);
                    });
                }
            } // mfi
        } // vdir
    }
}
#endif
