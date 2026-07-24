#include "ERF.H"
#include "ERF_Utils.H"
#include "ERF_TerrainPoisson.H"

#include <AMReX_GMRES.H>

using namespace amrex;

/**
 * Solve the Poisson equation using FFT-preconditioned GMRES
 */
void ERF::solve_with_gmres (int lev, const Box& subdomain, MultiFab& rhs, MultiFab& phi,
                            Array<MultiFab,AMREX_SPACEDIM>& fluxes,
                            MultiFab& ax_sub, MultiFab& ay_sub, MultiFab& az_sub,
                            MultiFab& dJ_sub, MultiFab& znd_sub)
{
#ifdef ERF_USE_FFT
    BL_PROFILE("ERF::solve_with_gmres()");

    Real reltol = solverChoice.poisson_reltol;
    Real abstol = solverChoice.poisson_abstol;

    auto const dom_lo = lbound(Geom(lev).Domain());
    auto const dom_hi = ubound(Geom(lev).Domain());

    auto const sub_lo = lbound(subdomain);
    auto const sub_hi = ubound(subdomain);

    auto dx    = Geom(lev).CellSizeArray();

    Geometry my_geom;

    Array<int,AMREX_SPACEDIM> is_per; is_per[0] = 0; is_per[1] = 0; is_per[2] = 0;
    if (Geom(lev).isPeriodic(0) && sub_lo.x == dom_lo.x && sub_hi.x == dom_hi.x) { is_per[0] = 1;}
    if (Geom(lev).isPeriodic(1) && sub_lo.y == dom_lo.y && sub_hi.y == dom_hi.y) { is_per[1] = 1;}

    int coord_sys = 0;

    // If subdomain == domain then we pass Geom(lev) to the FFT solver
    if (subdomain == Geom(lev).Domain()) {
        my_geom.define(Geom(lev).Domain(), Geom(lev).ProbDomain(), coord_sys, is_per);
    } else {
        // else we create a new geometry based only on the subdomain
        // The information in my_geom used by the FFT routines is:
        //   1) my_geom.Domain()
        //   2) my_geom.CellSize()
        //   3) my_geom.isAllPeriodic() / my_geom.periodicity()
        RealBox rb( sub_lo.x   *dx[0],  sub_lo.y   *dx[1],  sub_lo.z   *dx[2],
                   (sub_hi.x+1)*dx[0], (sub_hi.y+1)*dx[1], (sub_hi.z+1)*dx[2]);
        my_geom.define(subdomain, rb, coord_sys, is_per);
    }

    amrex::GMRES<MultiFab, TerrainPoisson> gmsolver;

    TerrainPoisson tp(my_geom, rhs.boxArray(), rhs.DistributionMap(), domain_bc_type,
                      stretched_dz_d[lev], ax_sub, ay_sub, az_sub, dJ_sub, &znd_sub,
                      solverChoice.use_real_bcs);

    gmsolver.define(tp);

    gmsolver.setVerbose(mg_verbose);

    gmsolver.setRestartLength(50);

    static const int consistency_test = [] {
        int b = 0; ParmParse pp("erf"); pp.query("poisson_consistency_test", b); return b;
    }();
    // Mode 4 (probe b): run WITHOUT the FFT preconditioner to test whether
    // the GMRES recurrence is honest for the bare operator.
    tp.usePrecond(consistency_test == 4 ? false : true);

    // erf.poisson_consistency_test = 1: synthetic-phi check that the
    // operator (TerrainPoisson::apply) and the applied-correction path
    // (getFluxes -> momenta increment -> compute_divergence) represent the
    // SAME discrete operator. They must agree to round-off; the residual
    // divergence left after "converged" anelastic solves says they do not.
    // Writes plotfile "plt_poisson_consistency" with [Aphi, divflux, diff]
    // and aborts.
    if (consistency_test == 1) {
        const Box& dom = Geom(lev).Domain();
        const auto dlo = lbound(dom);
        const auto dhi = ubound(dom);
        MultiFab phi_t(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        for (MFIter mfi(phi_t); mfi.isValid(); ++mfi) {
            const Box& gbx = mfi.growntilebox(1);
            const Array4<Real>& arr = phi_t.array(mfi);
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real X = Real(i-dlo.x) / Real(dhi.x-dlo.x+1);
                Real Y = Real(j-dlo.y) / Real(dhi.y-dlo.y+1);
                Real Z = Real(k-dlo.z) / Real(dhi.z-dlo.z+1);
                arr(i,j,k) = std::sin(Real(2.0)*PI*X) * std::cos(Real(2.0)*PI*Y)
                           * std::cos(PI*Z) * Real(1000.0);
            });
        }
        MultiFab Aphi(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        tp.apply(Aphi, phi_t);

        Array<MultiFab,AMREX_SPACEDIM> flx{
            MultiFab(convert(rhs.boxArray(), IntVect(1,0,0)), rhs.DistributionMap(), 1, 0),
            MultiFab(convert(rhs.boxArray(), IntVect(0,1,0)), rhs.DistributionMap(), 1, 0),
            MultiFab(convert(rhs.boxArray(), IntVect(0,0,1)), rhs.DistributionMap(), 1, 0)};
        tp.getFluxes(phi_t, flx);

        MultiFab divf(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        Array<MultiFab const*,AMREX_SPACEDIM> fp{&flx[0], &flx[1], &flx[2]};
        compute_divergence(lev, divf, fp, Geom(lev));

        MultiFab diff(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab::Copy(diff, divf, 0, 0, 1, 0);
        MultiFab::Subtract(diff, Aphi, 0, 0, 1, 0);

        Print() << "POISSON CONSISTENCY: |Aphi| max/L2 = " << Aphi.norm0() << " " << Aphi.norm2()
                << "  |divflux| max/L2 = " << divf.norm0() << " " << divf.norm2()
                << "  |diff| max/L2 = " << diff.norm0() << " " << diff.norm2() << std::endl;

        MultiFab out(rhs.boxArray(), rhs.DistributionMap(), 3, 0);
        MultiFab::Copy(out, Aphi, 0, 0, 1, 0);
        MultiFab::Copy(out, divf, 0, 1, 1, 0);
        MultiFab::Copy(out, diff, 0, 2, 1, 0);
        WriteSingleLevelPlotfile("plt_poisson_consistency", out,
                                 {"Aphi", "divflux", "diff"}, Geom(lev), Real(0.0), 0);
        Abort("poisson_consistency_test complete -- see plt_poisson_consistency");
    }

    // Mode 3: identify the discrete left-null vector of A (probe c).
    // For candidate weight vectors w in {1, dJ}: if <A x, w> = 0 for
    // arbitrary x, then w spans null(A^T) and solvability of A phi = rhs
    // requires <rhs, w> = 0. The production mean subtraction zeroes the
    // dJ-weighted sum; if the actual null vector is the constant vector,
    // the incompatible remainder is the stall floor.
    if (consistency_test == 3) {
        const Box& dom3 = Geom(lev).Domain();
        const auto dlo3 = lbound(dom3);
        const auto dhi3 = ubound(dom3);
        const Real Ncells = Real(dom3.numPts());
        for (int trial = 0; trial < 2; ++trial) {
            MultiFab phi_t(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
            for (MFIter mfi(phi_t); mfi.isValid(); ++mfi) {
                const Box& gbx = mfi.growntilebox(1);
                const Array4<Real>& arr = phi_t.array(mfi);
                const Real ph = (trial == 0) ? Real(0.0) : Real(1.234);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    Real X = Real(i-dlo3.x) / Real(dhi3.x-dlo3.x+1);
                    Real Y = Real(j-dlo3.y) / Real(dhi3.y-dlo3.y+1);
                    Real Z = Real(k-dlo3.z) / Real(dhi3.z-dlo3.z+1);
                    arr(i,j,k) = std::sin(Real(2.0)*PI*X + ph) * std::cos(Real(4.0)*PI*Y - ph)
                               * std::cos(PI*Z + myhalf*ph) * Real(1000.0);
                });
            }
            MultiFab Aphi(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
            tp.apply(Aphi, phi_t);
            Real s_std = Aphi.sum();
            Real s_dJ  = MultiFab::Dot(Aphi, 0, dJ_sub, 0, 1, 0);
            Real nA    = Aphi.norm2();
            Real ndJ   = dJ_sub.norm2();
            Print() << "NULLSPACE trial " << trial
                    << ": <Ax,1>/(|Ax| sqrtN) = " << s_std/(nA*std::sqrt(Ncells))
                    << "   <Ax,dJ>/(|Ax| |dJ|) = " << s_dJ/(nA*ndJ) << std::endl;
        }
        Real r_std = rhs.sum();
        Real r_dJ  = MultiFab::Dot(rhs, 0, dJ_sub, 0, 1, 0);
        Real nr    = rhs.norm2();
        Print() << "NULLSPACE rhs: <rhs,1>/(|rhs| sqrtN) = " << r_std/(nr*std::sqrt(Ncells))
                << "   <rhs,dJ>/(|rhs| |dJ|) = " << r_dJ/(nr*dJ_sub.norm2()) << std::endl;
    }

    // Mode 5: determinism + linearity of the preconditioner (and operator).
    // GMRES assumes M^-1 is a fixed linear operator; if precond(v) differs
    // between calls or fails superposition, the recurrence is fiction.
    if (consistency_test == 5) {
        MultiFab v1(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        MultiFab v2(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        const Box& dom5 = Geom(lev).Domain();
        const auto dlo5 = lbound(dom5); const auto dhi5 = ubound(dom5);
        for (int t = 0; t < 2; ++t) {
            MultiFab& vv = (t==0) ? v1 : v2;
            for (MFIter mfi(vv); mfi.isValid(); ++mfi) {
                const Box& gbx = mfi.growntilebox(1);
                const Array4<Real>& arr = vv.array(mfi);
                const Real ph = (t==0) ? Real(0.4) : Real(2.7);
                ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                    Real X = Real(i-dlo5.x)/Real(dhi5.x-dlo5.x+1);
                    Real Y = Real(j-dlo5.y)/Real(dhi5.y-dlo5.y+1);
                    Real Z = Real(k-dlo5.z)/Real(dhi5.z-dlo5.z+1);
                    arr(i,j,k) = std::sin(Real(6.0)*PI*X+ph)*std::sin(Real(2.0)*PI*Y-ph)
                               * std::cos(Real(3.0)*PI*Z*Z+ph) * Real(100.0);
                });
            }
        }
        MultiFab o1(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        MultiFab o2(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        MultiFab o3(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        // Determinism: M^-1 v1 twice
        tp.precond(o1, v1); tp.precond(o2, v1);
        MultiFab::Subtract(o2, o1, 0, 0, 1, 0);
        Print() << "PRECOND DETERMINISM: |M v1| L2 = " << o1.norm2()
                << "  |diff| L2 = " << o2.norm2() << std::endl;
        // Linearity: M^-1(v1 + 2 v2) vs M^-1 v1 + 2 M^-1 v2
        MultiFab vsum(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        MultiFab::LinComb(vsum, one, v1, 0, Real(2.0), v2, 0, 0, 1, 1);
        tp.precond(o3, vsum);
        tp.precond(o2, v2);
        MultiFab::Saxpy(o1, Real(2.0), o2, 0, 0, 1, 0);     // o1 = Mv1 + 2 Mv2
        MultiFab::Subtract(o3, o1, 0, 0, 1, 0);
        Print() << "PRECOND LINEARITY: |M(v1+2v2)| basis L2 = " << o1.norm2()
                << "  |mismatch| L2 = " << o3.norm2() << std::endl;
        // Operator determinism for completeness
        MultiFab a1(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab a2(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        tp.apply(a1, v1); tp.apply(a2, v1);
        MultiFab::Subtract(a2, a1, 0, 0, 1, 0);
        Print() << "OPERATOR DETERMINISM: |A v1| L2 = " << a1.norm2()
                << "  |diff| L2 = " << a2.norm2() << std::endl;
        Abort("poisson_consistency_test mode 5 complete");
    }

    // Prescribed-inflow projection rhs prep: with only the domain-boundary
    // faces masked, the left-null vector of the masked operator is full dJ
    // (interior fluxes telescope; masked faces contribute zero). The
    // production mean subtraction plus enforceInOutSolvability should
    // already have zeroed this component; remove any remainder exactly so
    // the deflated GMRES sees a compatible rhs.
    if (solverChoice.use_real_bcs) {
        Real c = MultiFab::Dot(rhs, 0, dJ_sub, 0, 1, 0) / MultiFab::Dot(dJ_sub, 0, dJ_sub, 0, 1, 0);
        MultiFab::Saxpy(rhs, -c, dJ_sub, 0, 0, 1, 0);
        if (mg_verbose > 0) {
            Print() << "Prescribed-inflow rhs prep: removed dJ-component c = "
                    << c << std::endl;
        }
    }

    gmsolver.solve(phi, rhs, reltol, abstol);

    // Iterative refinement on the TRUE residual. The GMRES recurrence
    // residual on this preconditioned terrain operator under-reports the
    // true residual by orders of magnitude (measured in double: reported
    // ~1e-8 relative, true 12% relative; the un-projected remainder then
    // pumps (rho theta) in anelastic runs at ~5 K/step near the lid).
    // Each re-solve on the true residual reduces it ~10-50x; iterate until
    // ||rhs - A phi|| actually meets the requested tolerance.
    {
        const int  max_refine = 10;
        const Real rhsnorm    = rhs.norm2();
        const Real target     = std::max(reltol * rhsnorm, abstol);
        MultiFab phi_g (rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        MultiFab Aphi_r(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab res_r (rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab dphi  (phi.boxArray(), phi.DistributionMap(), 1, phi.nGrowVect());
        Real truenorm = -one;
        int iref = 0;
        for (; iref < max_refine; ++iref) {
            MultiFab::Copy(phi_g, phi, 0, 0, 1, 0);
            tp.apply(Aphi_r, phi_g);
            MultiFab::Copy(res_r, rhs, 0, 0, 1, 0);
            MultiFab::Subtract(res_r, Aphi_r, 0, 0, 1, 0);
            truenorm = res_r.norm2();
            if (truenorm <= target) { break; }
            dphi.setVal(zero);
            gmsolver.solve(dphi, res_r, reltol, abstol);
            MultiFab::Add(phi, dphi, 0, 0, 1, 0);
            // Accept the pass only if it reduced the true residual
            // (measured: passes can diverge -- see UPSTREAM_ISSUES).
            MultiFab::Copy(phi_g, phi, 0, 0, 1, 0);
            tp.apply(Aphi_r, phi_g);
            MultiFab::Copy(res_r, rhs, 0, 0, 1, 0);
            MultiFab::Subtract(res_r, Aphi_r, 0, 0, 1, 0);
            if (res_r.norm2() >= truenorm) {
                MultiFab::Subtract(phi, dphi, 0, 0, 1, 0);   // revert
                ++iref;
                break;
            }
        }
        if (mg_verbose > 0) {
            Print() << "GMRES refinement: " << iref << " passes, true |rhs - A phi| L2 = "
                    << truenorm << " (target " << target << ")" << std::endl;
        }
        if (consistency_test == 3 && truenorm > target) {
            // Alignment of the stalled residual with candidate null vectors
            const Real Nc = Real(Geom(lev).Domain().numPts());
            Real c1 = res_r.sum() / (truenorm * std::sqrt(Nc));
            Real c2 = MultiFab::Dot(res_r, 0, dJ_sub, 0, 1, 0) / (truenorm * dJ_sub.norm2());
            Print() << "NULLSPACE stalled residual: cos(res,1) = " << c1
                    << "   cos(res,dJ) = " << c2 << std::endl;
        }
        if (truenorm > Real(100.0) * target) {
            Print() << "WARNING: GMRES refinement stalled: true residual " << truenorm
                    << " vs target " << target << std::endl;
        }
    }

    tp.getFluxes(phi, fluxes);

    // Mode 2: verify the SOLVED phi / real fluxes at the real point in the
    // sequence: || div(fluxes) + rhs || should match the GMRES residual.
    if (consistency_test == 2) {
        MultiFab divf(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        Array<MultiFab const*,AMREX_SPACEDIM> fp{&fluxes[0], &fluxes[1], &fluxes[2]};
        compute_divergence(lev, divf, fp, Geom(lev));
        MultiFab chk(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab::Copy(chk, divf, 0, 0, 1, 0);
        MultiFab::Add(chk, rhs, 0, 0, 1, 0);
        Print() << "SOLVED-PHI CONSISTENCY: |rhs| L2 = " << rhs.norm2()
                << "  |divflux| L2 = " << divf.norm2()
                << "  |divflux + rhs| max/L2 = " << chk.norm0() << " " << chk.norm2() << std::endl;

        // True operator residual: rhs - A phi_solved (discriminates GMRES
        // recurrence breakdown from flux-path linearity violation)
        MultiFab Aphi_s(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab phi_g(rhs.boxArray(), rhs.DistributionMap(), 1, 1);
        MultiFab::Copy(phi_g, phi, 0, 0, 1, 0);
        tp.apply(Aphi_s, phi_g);
        MultiFab res(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab::Copy(res, rhs, 0, 0, 1, 0);
        MultiFab::Subtract(res, Aphi_s, 0, 0, 1, 0);
        Print() << "TRUE OPERATOR RESIDUAL: |rhs - A phi| max/L2 = "
                << res.norm0() << " " << res.norm2() << std::endl;
    }

    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);
        const Array4<Real      >& fx_ar = fluxes[0].array(mfi);
        const Array4<Real      >& fy_ar = fluxes[1].array(mfi);
        const Array4<Real const>& mf_ux = mapfac[lev][MapFacType::u_x]->const_array(mfi);
        const Array4<Real const>& mf_vy = mapfac[lev][MapFacType::v_y]->const_array(mfi);
        ParallelFor(xbx,ybx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fx_ar(i,j,k) *= mf_ux(i,j,0);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fy_ar(i,j,k) *= mf_vy(i,j,0);
        });
    } // mfi
#else
    amrex::ignore_unused(lev, rhs, phi, fluxes, ax_sub, ay_sub, az_sub, dJ_sub, znd_sub);
#endif

    // ****************************************************************************
    // Impose bc's on pprime
    // ****************************************************************************
    ImposeBCsOnPhi(lev, phi, subdomain);
}
