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

    tp.usePrecond(true);

    // erf.poisson_consistency_test = 1: synthetic-phi check that the
    // operator (TerrainPoisson::apply) and the applied-correction path
    // (getFluxes -> momenta increment -> compute_divergence) represent the
    // SAME discrete operator. They must agree to round-off; the residual
    // divergence left after "converged" anelastic solves says they do not.
    // Writes plotfile "plt_poisson_consistency" with [Aphi, divflux, diff]
    // and aborts.
    static const int consistency_test = [] {
        int b = 0; ParmParse pp("erf"); pp.query("poisson_consistency_test", b); return b;
    }();
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
