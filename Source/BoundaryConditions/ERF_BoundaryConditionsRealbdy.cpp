#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_EOS.H"

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
    // DECOUPLED from hindcast_mass_consistent_bdy. That knob used to imply this one,
    // which made attempt 2 a compound experiment: band rho relaxation AND rho specified
    // at the wall face. The characteristic analysis says the wall setting is
    // inadmissible -- rho rides the OUTGOING u_n - c wave and must be computed from the
    // interior -- and attempt 2's failure was near-wall (excess at d<4 amplifying with
    // height), i.e. consistent with the wall half being the culprit. Band-only rho
    // relaxation with the wall left free is a different experiment and was never
    // separable until now. Only erf.hindcast_bdy_rho sets the wall.
    static const bool l_bdy_rho = [] {
        bool b=false; ParmParse pp("erf");
        pp.query("hindcast_bdy_rho", b);
        return b; }();

    // MPAS-A precedent: set w = 0 in the specified zone instead of taking a
    // zero-gradient copy from the interior. Reported to alleviate spurious
    // streamers and near-boundary instability under strong inflow.
    static const bool l_w_zero = [] {
        bool b=false; ParmParse pp("erf");
        pp.query("hindcast_bdy_w_zero", b); return b; }();
    const bool l_have_ext = (!bdy_data_xlo.empty()) &&
        (static_cast<int>(bdy_data_xlo[0].size()) > HindcastBdyVars::RHO);
    if (l_bdy_rho && l_have_ext) { cons_read[Rho_comp] = 1; }

    // Supply cloud and rain water at the lateral boundary instead of zeroing
    // them. The frame carries both and they were already being interpolated;
    // zeroing them forces the model to grow condensate from scratch over the
    // 120-165 km measured for this domain. NOTE: ice species (qi/qs/qg) are NOT
    // carried -- the frames hold no ciwc/cswc -- and Roberge et al. find ice
    // dominant for winter cases, so expect partial recovery, not full.
    static const bool l_bdy_hydro = [] {
        bool b=false; ParmParse pp("erf");
        pp.query("hindcast_bdy_hydrometeors", b); return b; }();
    const bool l_have_hydro = l_bdy_hydro && (!bdy_data_xlo.empty()) &&
        (static_cast<int>(bdy_data_xlo[0].size()) > HindcastBdyVars::QR);
    if (l_have_hydro) { cons_read[RhoQ2_comp] = 1; cons_read[RhoQ3_comp] = 1; }

    Vector<Vector<int>> is_read;
    is_read.push_back( cons_read );
    is_read.push_back( {1} ); // xvel
    is_read.push_back( {1} ); // yvel
    is_read.push_back( {0} ); // zvel

    // Real BC mapping (WRF/MetGrid)
    Vector<int> cons_map = {HindcastBdyVars::RHO, RealBdyVars::T, RhoKE_comp, RhoScalar_comp,
                            RealBdyVars::QV, HindcastBdyVars::QC, HindcastBdyVars::QR,
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
                    // Zero the moisture species we cannot supply. With
                    // hindcast_bdy_hydrometeors that is everything beyond rain
                    // (i.e. the ice species); without it, everything beyond qv.
                    const int last_supplied = l_have_hydro ? RhoQ3_comp : RhoQ1_comp;
                    const bool zero_here = (comp_idx > last_supplied) ||
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

    // ***********************************************************************************
    // CORNER w CONTAINMENT (erf.realbdy_corner_w_clamp, default on)
    //
    // The four cells where two lateral bands meet run away in |w|: measured 18.0 m/s
    // at (0,95) on 192x96 against 0.2-0.8 two cells along the same wall, growing
    // monotonically from 14 h and going non-finite at 18 h. On 128x64 the same four
    // cells hold the domain maximum (7.75 against an interior 2.09) on every run we
    // scored as clean. See UPSTREAM_ISSUES 27.
    //
    // This is CONTAINMENT, not a fix. Three candidate mechanisms were tested and
    // falsified, and the falsifications are what justify clamping rather than
    // continuing to hunt:
    //   * corner weight MAGNITUDE and the grad(F) branch discontinuity: replacing
    //     max(xi^2,eta^2) with a C1 blend changed F by up to 27% at individual corner
    //     cells and moved the response 0.6%, dying at 18.001 h against 18.002 h
    //     (item 27d). Near-zero sensitivity to F rules out the weighting, and with it
    //     min(F_x,F_y) and a corner taper, which act on the same coefficient.
    //   * double-writing / later-launch-wins: realbdy_bc_bxs_xy explicitly trims the
    //     y-face boxes in x ("Remove overlapping corners from y-face boxes"), so a
    //     corner is written exactly once, by the x box.
    //   * disagreeing per-face targets: all four planes are filled by
    //     strip_to_global_fab from the SAME MultiFab, so bdy_data_xlo and
    //     bdy_data_yhi at a corner are copies of one source element -- bit-identical,
    //     for every variable.
    //
    // The four corner cells are inside the discard region at any resolution (the outer
    // ~10 cells are forcing-dominated and excluded from every diagnostic), so clamping
    // |w| there costs nothing scientific. Clamp to the largest |w| among the 8
    // horizontal neighbours, preserving sign, so the corner can never lead the domain.
    if (!cons_only && corner_w_clamp()) {
        MultiFab& wmf = *mfs[Vars::zvel];
        const auto& dlo = lbound(geom[lev].Domain());
        const auto& dhi = ubound(geom[lev].Domain());
        const int ci[4] = {dlo.x, dlo.x, dhi.x, dhi.x};
        const int cj[4] = {dlo.y, dhi.y, dhi.y, dlo.y};
        for (MFIter mfi(wmf); mfi.isValid(); ++mfi) {
            const Box& vbx = mfi.validbox();
            const Array4<Real>& w = wmf.array(mfi);
            for (int c = 0; c < 4; ++c) {
                const int I = ci[c], J = cj[c];
                // Inward neighbour offsets: a corner has neighbours only on the
                // inward side, so step toward the interior in each direction.
                const int si = (I == dlo.x) ? 1 : -1;
                const int sj = (J == dlo.y) ? 1 : -1;
                Box col(IntVect(I,J,vbx.smallEnd(2)), IntVect(I,J,vbx.bigEnd(2)));
                col &= vbx;
                if (!col.ok()) { continue; }
                ParallelFor(col, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real nmax = Real(0.0);
                    for (int dj = 0; dj <= 2; ++dj) {
                        for (int di = 0; di <= 2; ++di) {
                            if (di == 0 && dj == 0) { continue; }
                            nmax = amrex::max(nmax, std::abs(w(i+si*di, j+sj*dj, k)));
                        }
                    }
                    const Real a = std::abs(w(i,j,k));
                    if (a > nmax) {
                        w(i,j,k) = (w(i,j,k) > Real(0.0)) ? nmax : -nmax;
                    }
                });
            }
        }
    }

    // ***********************************************************************************
    // NSCBC lateral boundary treatment (erf.nscbc_lateral)
    //
    // Replaces the "specify everything, everywhere" fill above with the admissible
    // characteristic set, determined per face point per level per BC application.
    //
    // Counting, for ERF's boundary-normal system (n = 5 + N variables, N scalars):
    //   lambda = u_n - c  (x1),  u_n  (x3+N),  u_n + c  (x1)
    // Subsonic INFLOW  : 4 + N admissible conditions -- u_n, both tangential
    //                    velocities, theta, and every scalar.  rho is left FREE:
    //                    it rides the OUTGOING u_n - c wave and must be computed
    //                    from the interior.
    // Subsonic OUTFLOW : exactly ONE, independent of N.
    //
    // The incoming acoustic characteristic, per face.  Working in the OUTWARD-normal
    // frame (u_n = u.nhat, nhat pointing out of the domain), the incoming acoustic
    // wave is always lambda = u_n - c < 0, carrying J- = u_n - 2c/(gamma-1).
    // Translated back to grid components that is:
    //
    //   face   nhat    u_n        outflow means   incoming acoustic   invariant specified
    //   xlo    -x      -u         u < 0           lambda_x = u + c    J+ = u + 2c/(g-1)
    //   xhi    +x      +u         u > 0           lambda_x = u - c    J- = u - 2c/(g-1)
    //   ylo    -y      -v         v < 0           lambda_y = v + c    J+ = v + 2c/(g-1)
    //   yhi    +y      +v         v > 0           lambda_y = v - c    J- = v - 2c/(g-1)
    //
    // i.e. the LO faces specify J+ and the HI faces specify J-.  Getting this backwards
    // would impose an OUTGOING invariant -- wrong, but it would still run.  The code
    // below works in the outward-normal frame so there is a single formula and the
    // per-face sign lives only in nsign.
    //
    // At INFLOW the outgoing wave is the other one (lambda_x = u - c at xlo), which is
    // what carries rho/p out of the domain -- hence rho free.
    //
    // APPROXIMATION 1 (isentropic invariants).  J+- = u_n +- 2c/(gamma-1) is the 1-D
    // ISENTROPIC form.  In a stratified atmosphere entropy varies along the normal so
    // these are not exactly conserved.  Over one acoustic substep (dtau ~ 0.4 s, sound
    // travels ~134 m against dx = 3 km) the error is small, and theta is carried by the
    // lambda_0 family which we take from the interior separately -- the isentropic
    // assumption enters only in the acoustic pair.
    //
    // APPROXIMATION 2 (the sigma blend).  Blending a (4+N)-condition inflow state with
    // a 1-condition outflow state is NOT a preserved count inside the blend window; it
    // is a Robin-type condition there. Marchesiello et al. instead blend the relaxation
    // TIMESCALE between inflow and outflow values, keeping the condition set fixed.
    // The blend is used here because it is pointwise in u_n rather than in distance from
    // the wall, so it cannot rebuild the monotone wall-normal ramp whose gradient is
    // the artifact -- but it is an approximation, not a derivation.  eps is small
    // (0.5 m/s) so sigma is 0 or 1 almost everywhere.  If it flaps, replace with a hard
    // switch plus hysteresis on a persistent per-face flag.
    // ***********************************************************************************
    static const int l_nscbc_outflow = [] {
        int v=0; ParmParse pp("erf"); pp.query("nscbc_outflow", v); return v; }();
    static const Real l_nscbc_eps = [] {
        Real e=Real(0.5); ParmParse pp("erf"); pp.query("nscbc_eps", e); return e; }();
    // Bisection bitmask: 1=cons pass, 2=velocity pass, 4=specify KE/scalar at
    // inflow, 8=specify w=0 at inflow, 16=characteristic inflow density.
    // Default 31 = the full formulation.
    static const int l_nscbc_parts = [] {
        int v=31; ParmParse pp("erf"); pp.query("nscbc_parts", v); return v; }();
    static const Real l_nscbc_tke = [] {
        Real t=Real(0.01); ParmParse pp("erf"); pp.query("nscbc_inflow_tke", t); return t; }();

    if (nscbc_lateral() && !cons_only)
    {
        MultiFab& cons_mf = *mfs[Vars::cons];
        MultiFab& xvel_mf = *mfs[Vars::xvel];
        MultiFab& yvel_mf = *mfs[Vars::yvel];
        MultiFab& zvel_mf = *mfs[Vars::zvel];

        // `domain` and `set_width` above are scoped to the per-variable loop
        const Box domain_n   = geom[lev].Domain();
        const int set_width_n = 1;
        const auto dlo = lbound(domain_n);
        const auto dhi = ubound(domain_n);

        // Scratch carrying the per-column boundary solution, stored at the wall-adjacent
        // CELL so both the cons pass and the velocity passes read the same numbers:
        //   0 = sigma (1 = inflow, 0 = outflow)
        //   1 = rho at the boundary       (outflow, variant 1; else the extrapolated rho)
        //   2 = rho*theta at the boundary (   "   )
        //   3 = wall-normal velocity, OUTWARD-positive
        //   4 = inflow density ratio rho_characteristic / rho_zero-gradient
        const int ng_nsc = std::max(ngvect_cons.max(), ngvect_vels.max());
        MultiFab nsc(cons_mf.boxArray(), cons_mf.DistributionMap(), 5, ng_nsc);
        nsc.setVal(0.0);

        const Real gm1  = Gamma - Real(1.0);
        const int  ncmp = ncomp_cons;
        const int  icmp = icomp_cons;
        const int  ovar = l_nscbc_outflow;
        const Real eps  = l_nscbc_eps;
        const int  prt  = l_nscbc_parts;
        const Real tke_in = l_nscbc_tke;

        MultiFab r_hse(base_state[lev], make_alias, BaseState::r0_comp, 1);
        MultiFab p_hse(base_state[lev], make_alias, BaseState::p0_comp, 1);

        // fdir = 0:xlo 1:xhi 2:ylo 3:yhi
        for (int fdir = 0; fdir < 4; ++fdir)
        {
            const int  ndir  = (fdir < 2) ? 0 : 1;                 // normal coordinate
            const Real nsign = (fdir % 2 == 0) ? Real(-1.) : Real(1.); // nhat . e_ndir
            // Wall-adjacent cell, first interior cell, and the wall / first-interior face
            const int  wcell = (fdir==0) ? dlo.x : (fdir==1) ? dhi.x : (fdir==2) ? dlo.y : dhi.y;
            const int  icell = wcell - static_cast<int>(nsign);    // one cell inward
            const int  wface = (nsign < 0) ? wcell : wcell + 1;    // the wall face index
            const int  iface = wface - static_cast<int>(nsign);    // first face inward

            // ---- Stage 1: solve for the boundary state, from the UNMODIFIED interior ----
            for (MFIter mfi(nsc,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Wall-adjacent slab, intersected with this tile (setSmall/setBig alone
                // would GROW a tile that does not touch the wall)
                Box wall_slab = domain_n;
                if (ndir == 0) { wall_slab.setSmall(0,wcell); wall_slab.setBig(0,wcell); }
                else           { wall_slab.setSmall(1,wcell); wall_slab.setBig(1,wcell); }
                Box bx = mfi.tilebox() & wall_slab;
                if (bx.isEmpty()) continue;

                const Array4<Real>&       ns = nsc.array(mfi);
                const Array4<const Real>& cs = cons_mf.const_array(mfi);
                const Array4<const Real>& uu = xvel_mf.const_array(mfi);
                const Array4<const Real>& vv = yvel_mf.const_array(mfi);
                const Array4<const Real>& rh = r_hse.const_array(mfi);
                const Array4<const Real>& ph = p_hse.const_array(mfi);

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    // Driver normal velocity at the wall face, made OUTWARD-positive.
                    // This is the just-prescribed boundary value (the Eta rule): it is
                    // smooth in time, so the regime does not chatter on model noise.
                    Real un_w = (ndir == 0) ? nsign * uu(wface,j,k) : nsign * vv(i,wface,k);

                    // sigma = 1 deep inflow (un_w <= -eps), 0 deep outflow (un_w >= +eps)
                    Real sig = Real(0.5) - un_w / (Real(2.0)*eps);
                    sig = amrex::min(amrex::max(sig, Real(0.0)), Real(1.0));
                    ns(i,j,k,0) = sig;

                    // Interior state, one cell / one face in
                    const int ii = (ndir == 0) ? icell : i;
                    const int jj = (ndir == 1) ? icell : j;
                    Real rho_i = cs(ii,jj,k,Rho_comp);
                    Real rt_i  = cs(ii,jj,k,RhoTheta_comp);
                    Real un_i  = (ndir == 0) ? nsign * uu(iface,j,k) : nsign * vv(i,iface,k);

                    // Default (variant 0): zeroth-order extrapolation from the interior,
                    // which for an outflow face IS the upwind-biased stencil.
                    Real rho_b = rho_i, rt_b = rt_i, un_b = un_i;

                    if (ovar == 1) {
                        Real p_i = getPgivenRTh(rt_i);
                        Real c_i = std::sqrt(Gamma * p_i / rho_i);
                        // Only meaningful subsonic; fall back to extrapolation otherwise.
                        if (c_i > Real(0.0) && std::abs(un_i) < c_i) {
                            // Outgoing invariant, carried out of the domain by lambda = u_n + c
                            Real Jout = un_i + Real(2.0)*c_i/gm1;
                            // The ONE incoming condition: J- from the far field.  p_inf and
                            // rho_inf are the hydrostatic base state; u_n,inf is the driver.
                            Real c_inf = std::sqrt(Gamma * ph(i,j,k) / rh(i,j,k));
                            Real Jin   = un_w - Real(2.0)*c_inf/gm1;

                            Real un_s = Real(0.5)*(Jout + Jin);
                            Real c_s  = Real(0.25)*gm1*(Jout - Jin);
                            if (c_s > Real(0.0)) {
                                // Entropy comes from the interior (outgoing lambda_0 wave)
                                Real Kent = p_i / std::pow(rho_i, Gamma);
                                Real rho_s = std::pow(c_s*c_s/(Gamma*Kent), Real(1.0)/gm1);
                                Real p_s   = Kent * std::pow(rho_s, Gamma);
                                rho_b = rho_s;
                                rt_b  = getRhoThetagivenP(p_s);
                                un_b  = un_s;
                            }
                        }
                    }
                    ns(i,j,k,1) = rho_b;
                    ns(i,j,k,2) = rt_b;
                    ns(i,j,k,3) = un_b;

                    // INFLOW density, from the OUTGOING acoustic characteristic.
                    // In the outward-normal frame u_n + c is outgoing at a subsonic
                    // lateral boundary in BOTH regimes, so J+ = u_n + 2c/(gamma-1)
                    // carries the interior state out and fixes the boundary sound
                    // speed once u_n is specified:  c_b = (gamma-1)/2 * (J+ - u_n,drv).
                    // With theta also specified, rho follows from
                    //   c^2 = gamma * theta * p_0 * (R_d/p_0)^gamma * (rho theta)^(gamma-1)
                    // taken as a RATIO against the interior, which needs no constants:
                    //   rho_b/rho_i = (th_i/th_b) * [ (c_b/c_i)^2 (th_i/th_b) ]^(1/(gamma-1))
                    // Zero-gradient rho -- what ERF does today -- is the crudest possible
                    // stand-in for this, and it puts a density discontinuity on the wall
                    // face, which is exactly where the global flux sum is evaluated.
                    Real ratio = Real(1.0);
                    Real p_ii  = getPgivenRTh(rt_i);
                    Real c_ii  = std::sqrt(Gamma * p_ii / rho_i);
                    Real th_i  = rt_i / rho_i;
                    Real rho_w = cs(i,j,k,Rho_comp);
                    Real th_b  = (rho_w > Real(0.0)) ? cs(i,j,k,RhoTheta_comp)/rho_w : th_i;
                    Real Jp    = un_i + Real(2.0)*c_ii/gm1;
                    Real c_bb  = Real(0.5)*gm1*(Jp - un_w);
                    if (c_bb > Real(0.0) && c_ii > Real(0.0) && th_b > Real(0.0)) {
                        Real r  = th_i / th_b;
                        Real cr = (c_bb/c_ii)*(c_bb/c_ii);
                        ratio = r * std::pow(cr*r, Real(1.0)/gm1);
                    }
                    if (!(prt & 16)) ratio = Real(1.0);
                    ns(i,j,k,4) = amrex::min(amrex::max(ratio, Real(0.5)), Real(2.0));
                });
            } // mfi
            nsc.FillBoundary(geom[lev].periodicity());

            // ---- Stage 2: apply to cons in the specified zone + exterior ghosts ----
            for (MFIter mfi(cons_mf,TilingIfNotGPU()); (prt & 1) && mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(ngvect_cons);
                Box b_xlo, b_xhi, b_ylo, b_yhi;
                realbdy_bc_bxs_xy(gbx, domain_n, set_width_n, b_xlo, b_xhi, b_ylo, b_yhi, ngvect_cons);
                Box bx = (fdir==0) ? b_xlo : (fdir==1) ? b_xhi : (fdir==2) ? b_ylo : b_yhi;
                if (bx.isEmpty()) continue;

                const Array4<Real>&       cs = cons_mf.array(mfi);
                const Array4<const Real>& ns = nsc.const_array(mfi);

                // Rho FIRST, in its own kernel: the scalar conditions below need the
                // updated boundary density, and reading it in the same kernel that
                // writes it would be a race.
                if (icmp == 0) {
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const int wi = (ndir == 0) ? wcell : amrex::min(amrex::max(i,dlo.x),dhi.x);
                        const int wj = (ndir == 1) ? wcell : amrex::min(amrex::max(j,dlo.y),dhi.y);
                        Real sig = ns(wi,wj,k,0);
                        // rho is FREE at inflow -- the zero-gradient value already in place
                        // IS the statement that it is computed from the interior.
                        cs(i,j,k,Rho_comp) = sig*cs(i,j,k,Rho_comp)*ns(wi,wj,k,4)
                                           + (Real(1.0)-sig)*ns(wi,wj,k,1);
                    });
                }

                ParallelFor(bx, ncmp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    const int c = n + icmp;
                    if (c == Rho_comp) return;
                    // Index of the wall-adjacent cell for this column
                    const int wi = (ndir == 0) ? wcell : amrex::min(amrex::max(i,dlo.x),dhi.x);
                    const int wj = (ndir == 1) ? wcell : amrex::min(amrex::max(j,dlo.y),dhi.y);
                    const int ii = (ndir == 0) ? icell : wi;
                    const int jj = (ndir == 1) ? icell : wj;

                    Real sig = ns(wi,wj,k,0);

                    // INFLOW value.  rho must stay FREE -- the fill above already left it
                    // zero-gradient, which is the correct (zeroth-order) statement that it
                    // is computed from the interior.  KE and the passive scalar are
                    // admissible inflow conditions that the fill above left floating;
                    // specify them (laminar, tracer-free inflow).  Everything else --
                    // theta, q_v, and the zeroed hydrometeors -- is already specified.
                    // rho-weighted quantities were built on the zero-gradient rho by the
                    // driver fill; rescale them onto the characteristic rho so theta and
                    // every q keep their SPECIFIED values.
                    Real v_in = cs(i,j,k,c) * ns(wi,wj,k,4);
                    if (prt & 4) {
                        // KE and the passive scalar are admissible inflow conditions that
                        // the driver fill leaves floating.  Specify them: a tracer-free
                        // inflow, and free-stream TKE at a small positive floor (a hard
                        // zero is pathological for MYNN's length-scale closure).
                        if      (c == RhoKE_comp)     v_in = cs(i,j,k,Rho_comp) * tke_in;
                        else if (c == RhoScalar_comp) v_in = Real(0.0);
                    }

                    // OUTFLOW value: computed from the interior.  rho and rho*theta take
                    // the characteristic solution when variant 1 is active.
                    Real v_out;
                    if      (c == Rho_comp)      { v_out = ns(wi,wj,k,1); }
                    else if (c == RhoTheta_comp) { v_out = ns(wi,wj,k,2); }
                    else                         { v_out = cs(ii,jj,k,c); }

                    cs(i,j,k,c) = sig*v_in + (Real(1.0)-sig)*v_out;
                });
            } // mfi

            // ---- Stage 3: apply to the velocities ----
            for (int vd = 0; vd < 3; ++vd)
            {
                if (!(prt & 2)) continue;
                MultiFab& vmf = (vd==0) ? xvel_mf : (vd==1) ? yvel_mf : zvel_mf;
                Box domv = domain_n; domv.convert(vmf.boxArray().ixType());

                for (MFIter mfi(vmf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    Box gbx = mfi.growntilebox(ngvect_vels);
                    Box b_xlo, b_xhi, b_ylo, b_yhi;
                    realbdy_bc_bxs_xy(gbx, domv, set_width_n, b_xlo, b_xhi, b_ylo, b_yhi, ngvect_vels);
                    Box bx = (fdir==0) ? b_xlo : (fdir==1) ? b_xhi : (fdir==2) ? b_ylo : b_yhi;
                    if (bx.isEmpty()) continue;

                    const Array4<Real>&       vl = vmf.array(mfi);
                    const Array4<const Real>& ns = nsc.const_array(mfi);
                    const bool is_normal = (vd == ndir);

                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const int wi = (ndir == 0) ? wcell : amrex::min(amrex::max(i,dlo.x),dhi.x);
                        const int wj = (ndir == 1) ? wcell : amrex::min(amrex::max(j,dlo.y),dhi.y);
                        Real sig = ns(wi,wj,k,0);

                        // INFLOW: u_n and the tangential horizontal velocity are already the
                        // driver values.  w is a boundary-TANGENTIAL velocity here and is an
                        // admissible inflow condition, but the driver carries no usable w
                        // (the frame slot is ERA5 omega in Pa/s), so specify w = 0.  The
                        // ERA5-consistent value is ~6 mm/s (measured), far below anything
                        // this run resolves.
                        Real v_in = ((vd == 2) && (prt & 8)) ? Real(0.0) : vl(i,j,k);

                        // OUTFLOW: computed from the interior.  The normal component takes
                        // the characteristic solution; tangential and w extrapolate.
                        Real v_out;
                        if (is_normal) {
                            v_out = nsign * ns(wi,wj,k,3);
                        } else {
                            const int ii = (ndir == 0) ? icell : i;
                            const int jj = (ndir == 1) ? icell : j;
                            v_out = vl(ii,jj,k);
                        }
                        vl(i,j,k) = sig*v_in + (Real(1.0)-sig)*v_out;
                    });
                } // mfi
            } // vd
        } // fdir
    }

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

    // ***********************************************************************************
    // Global mass constraint for the NSCBC lateral path (erf.nscbc_mass_tau)
    //
    // The characteristic count says how many conditions are ADMISSIBLE; it says nothing
    // about the discrete global mass budget.  In a limited-area domain that has to be
    // imposed separately, which is standard practice rather than invention: Flather
    // (1976) is derived from mass conservation, and ROMS/NEMO pair radiation OBCs with
    // an explicit barotropic inflow-outflow adjustment to preserve total volume.  ERF
    // already does the same thing for the anelastic path in enforceInOutSolvability.
    //
    // Here the total wall mass flux is rescaled to the value the ERA5 column-mass
    // tendency asks for.  That target is trustworthy: the ERA5 field is mass-consistent
    // under ERF's OWN discrete operator to an equivalent |w| of 6.4 mm/s (measured by
    // hindcast_check_mass_consistency), so we are rescaling toward a known small number.
    //
    //   Phi_desired = (M - M_tgt) / tau        (net OUTWARD flux; M > M_tgt => expel)
    //   du          = (Phi_desired - Phi_now) / sum_outflow(rho*A)
    //   u_n -> u_n + du        on OUTFLOW faces only
    //
    // The correction is ADDITIVE, not multiplicative.  A multiplicative rescaling was
    // tried first and is unusable: with Phi_desired ~ 0 it demands lambda = -Phi_in/Phi_out,
    // which in a net-convergent synoptic regime is many-fold, and multiplying outflow
    // velocities by that factor collapsed the CFL and killed the run inside 30 steps.
    // The additive form is what ROMS/NEMO actually apply to the barotropic mode, and its
    // magnitude is set by the imbalance divided by the wall's mass-flux capacity -- of
    // order 0.6 m/s for the +96 %/day drift being corrected here.
    //
    // du is a single scalar applied uniformly.  It deliberately has NO wall-normal
    // structure -- a spatially varying wall-normal correction would rebuild exactly the
    // grad(F) gradient this scheme exists to remove.  With set_width = 1 there is no
    // wall-normal extent to ramp over in any case.
    //
    // Runs last so it has the final word on the wall-normal velocity, after both the
    // characteristic pass and the barotropic wall flux correction.
    // ***********************************************************************************
    static const Real l_mass_tau = [] {
        Real t=Real(-1.0); ParmParse pp("erf"); pp.query("nscbc_mass_tau", t); return t; }();

    if (nscbc_lateral() && (l_mass_tau > Real(0.0)) && !cons_only &&
        solverChoice.init_type == InitType::HindCast &&
        !forecast_state_interp[lev].empty())
    {
        MultiFab& cons_mf = *mfs[Vars::cons];
        const MultiFab& tgt_mf = forecast_state_interp[lev][Vars::cons];
        const Box  dom_g = geom[lev].Domain();
        const auto glo = lbound(dom_g);
        const auto ghi = ubound(dom_g);
        const Real dxg = geom[lev].CellSize(0);
        const Real dyg = geom[lev].CellSize(1);
        const Real dzg = geom[lev].CellSize(2);

        // Total model and target mass.  Volume = detJ * dx*dy*dz (terrain-following
        // Jacobian); map factors are ~1.00 over this domain and are neglected, as they
        // are in compute_wall_flux_correction.
        const Real cellvol = dxg*dyg*dzg;
        Real M    = cellvol * amrex::ReduceSum(cons_mf, *detJ_cc[lev], 0,
                        [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                                   Array4<const Real> const& c,
                                                   Array4<const Real> const& J) -> Real {
                            Real s = 0.0;
                            AMREX_LOOP_3D(bx, i, j, k, { s += c(i,j,k,Rho_comp)*J(i,j,k); });
                            return s; });
        Real Mt   = cellvol * amrex::ReduceSum(tgt_mf, *detJ_cc[lev], 0,
                        [=] AMREX_GPU_HOST_DEVICE (Box const& bx,
                                                   Array4<const Real> const& c,
                                                   Array4<const Real> const& J) -> Real {
                            Real s = 0.0;
                            AMREX_LOOP_3D(bx, i, j, k, { s += c(i,j,k,Rho_comp)*J(i,j,k); });
                            return s; });

        // Wall mass flux, split into inward and outward parts.  Face area for an x-wall
        // is ax*dy*dz; corner cells own two walls and contribute to both.
        MultiFab wf(cons_mf.boxArray(), cons_mf.DistributionMap(), 3, 0);
        // (cons ghosts were filled by the driver/NSCBC passes above)
        wf.setVal(0.0);
        for (MFIter mfi(wf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Array4<Real>&       w  = wf.array(mfi);
            const Array4<const Real>& cs = cons_mf.const_array(mfi);
            const Array4<const Real>& uu = (*mfs[Vars::xvel]).const_array(mfi);
            const Array4<const Real>& vv = (*mfs[Vars::yvel]).const_array(mfi);
            const Array4<const Real>& axa = ax[lev]->const_array(mfi);
            const Array4<const Real>& aya = ay[lev]->const_array(mfi);

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Face density, formed the SAME way the dycore forms it in
                // VelocityToMomentum -- the average of the two adjacent cells, not the
                // wall cell alone.  ax/ay already carry the vertical stretching
                // (ax = 0.5*(z_nd(k+1)-z_nd(k))/dz_ref), so ax*dy*dz_ref is the true
                // face area and no separate dz correction is needed.
                Real fsum = Real(0.0), csum = Real(0.0);
                Real rf, f, a;
                if (i == glo.x) { rf = Real(0.5)*(cs(glo.x,j,k,Rho_comp)+cs(glo.x-1,j,k,Rho_comp));
                                  a = axa(glo.x  ,j,k)*dyg*dzg; f = rf*(-uu(glo.x  ,j,k))*a;
                                  fsum += f; if (f > Real(0.0)) csum += rf*a; }
                if (i == ghi.x) { rf = Real(0.5)*(cs(ghi.x,j,k,Rho_comp)+cs(ghi.x+1,j,k,Rho_comp));
                                  a = axa(ghi.x+1,j,k)*dyg*dzg; f = rf*( uu(ghi.x+1,j,k))*a;
                                  fsum += f; if (f > Real(0.0)) csum += rf*a; }
                if (j == glo.y) { rf = Real(0.5)*(cs(i,glo.y,k,Rho_comp)+cs(i,glo.y-1,k,Rho_comp));
                                  a = aya(i,glo.y  ,k)*dxg*dzg; f = rf*(-vv(i,glo.y  ,k))*a;
                                  fsum += f; if (f > Real(0.0)) csum += rf*a; }
                if (j == ghi.y) { rf = Real(0.5)*(cs(i,ghi.y,k,Rho_comp)+cs(i,ghi.y+1,k,Rho_comp));
                                  a = aya(i,ghi.y+1,k)*dxg*dzg; f = rf*( vv(i,ghi.y+1,k))*a;
                                  fsum += f; if (f > Real(0.0)) csum += rf*a; }
                w(i,j,k,0) = amrex::min(fsum, Real(0.0));   // inward  (negative)
                w(i,j,k,1) = amrex::max(fsum, Real(0.0));   // outward (positive)
                w(i,j,k,2) = csum;                          // rho*A on outflow faces
            });
        }
        Real Phi_in  = wf.sum(0);
        Real Phi_out = wf.sum(1);

        // Mass-flux capacity of the outflow faces: sum of rho*A over faces currently
        // carrying mass out.  This converts a flux deficit into a velocity increment.
        Real cap = wf.sum(2);
        Real Phi_des = (M - Mt) / l_mass_tau;
        Real du = Real(0.0);
        if (cap > Real(0.0)) { du = (Phi_des - (Phi_in + Phi_out)) / cap; }
        du = amrex::min(amrex::max(du, Real(-2.0)), Real(2.0));

        for (int vdir = 0; vdir < 2; ++vdir)
        {
            MultiFab& mf = (vdir == 0) ? *mfs[Vars::xvel] : *mfs[Vars::yvel];
            Box domv = dom_g; domv.convert(mf.boxArray().ixType());
            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                Box gbx = mfi.growntilebox(ngvect_vels);
                Box b_xlo, b_xhi, b_ylo, b_yhi;
                realbdy_bc_bxs_xy(gbx, domv, 1, b_xlo, b_xhi, b_ylo, b_yhi, ngvect_vels);
                const Array4<Real>& vel = mf.array(mfi);
                // Scale OUTFLOW faces only, by a single uniform factor.
                if (vdir == 0) {
                    ParallelFor(b_xlo, b_xhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { if (vel(i,j,k) < Real(0.0)) vel(i,j,k) -= du; },   // outward at xlo is -x
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { if (vel(i,j,k) > Real(0.0)) vel(i,j,k) += du; });
                } else {
                    ParallelFor(b_ylo, b_yhi,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { if (vel(i,j,k) < Real(0.0)) vel(i,j,k) -= du; },
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { if (vel(i,j,k) > Real(0.0)) vel(i,j,k) += du; });
                }
            }
        }
    }
}
#endif
