#include <ERF.H>
#include <ERF_Utils.H>

#ifdef ERF_USE_WINDFARM
#include <ERF_WindFarm.H>
#endif

using namespace amrex;

// Debug probe: erf.debug_cell = i j k prints rho, rhotheta, rhoQ1, rhoQ2 at
// that cell at each stage of ERF::Advance (budget attribution for the
// hindcast qv-pump investigation).
static void debug_cell_print (const MultiFab& S, const char* tag)
{
    static IntVect dbg(std::numeric_limits<int>::lowest(),0,0);
    static bool parsed = false;
    if (!parsed) {
        parsed = true;
        ParmParse pp("erf");
        Vector<int> v;
        if (pp.queryarr("debug_cell", v) && v.size() == 3) { dbg = IntVect(v[0],v[1],v[2]); }
    }
    if (dbg[0] == std::numeric_limits<int>::lowest()) { return; }
    Vector<Real> vals = get_cell_data(S, dbg);
    if (!vals.empty()) {
        const Real r = vals[Rho_comp];
        std::cout << "DBGCELL " << tag
                  << " rho "  << r
                  << " th "   << vals[RhoTheta_comp]/r
                  << " qv "   << ((static_cast<int>(vals.size()) > RhoQ1_comp) ? vals[RhoQ1_comp]/r : Real(0.))
                  << " qc "   << ((static_cast<int>(vals.size()) > RhoQ2_comp) ? vals[RhoQ2_comp]/r : Real(0.))
                  << std::endl;
    }
}

// Bound the prognostic MYNN TKE like WRF does. The MYNN-EDMF column solver
// clips its internal qke to <= 150 m2/s2 (2*TKE), and in WRF that clipped
// qke IS the prognostic state; in this port the dycore-advected RhoKE never
// sees the clip, and upper-level shear production (2.5-km cells at 8-15 km)
// was measured running TKE past 200 m2/s2 to NaN at ~4.6 h.
// (Free function: CUDA extended lambdas are not allowed inside the private
// member function ERF::Advance.)
// WRF-style specified-zone treatment for MASS in the hindcast lateral bands:
// blend density directly toward the interpolated forecast after each step,
// preserving theta and the moisture mixing ratio (the rho-weighted components
// are rescaled). Rationale (all measured at 3 km): the momentum-only RHS
// sponge accumulates mass along the whole band (~0.3 kg/m3 per hour; worst at
// band-overlap corners) until EOS/radiation blowup at ~4.6 h; an RHS
// mass-equation relaxation destabilizes the acoustic substepping in ~60
// steps. A direct post-step state blend does not enter the acoustic solver.
static void hindcast_blend_band_density (MultiFab& S, const MultiFab& fcons,
                                         const Geometry& geom, const Real dt,
                                         const SolverChoice& sc, const bool has_moisture)
{
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();
    const Real len      = sc.hindcast_lateral_sponge_length;
    const Real strength = sc.hindcast_lateral_sponge_strength;
    const Real xlo_end   = ProbLoArr[0] + len;
    const Real xhi_start = ProbHiArr[0] - len;
    const Real ylo_end   = ProbLoArr[1] + len;
    const Real yhi_start = ProbHiArr[1] - len;

    for (MFIter mfi(S); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        const Array4<Real>& s_arr = S.array(mfi);
        const Array4<Real const>& f_arr = fcons.const_array(mfi);
        const bool blend_q = has_moisture && (S.nComp() > RhoQ1_comp);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = ProbLoArr[0] + (i+myhalf) * dx[0];
            Real y = ProbLoArr[1] + (j+myhalf) * dx[1];
            Real xi = zero;
            if (x < xlo_end)   { xi = amrex::max(xi, (xlo_end - x)   / len); }
            if (x > xhi_start) { xi = amrex::max(xi, (x - xhi_start) / len); }
            if (y < ylo_end)   { xi = amrex::max(xi, (ylo_end - y)   / len); }
            if (y > yhi_start) { xi = amrex::max(xi, (y - yhi_start) / len); }
            if (xi > zero) {
                // Gentle rate: the measured pile-up is only ~1e-4 kg/m3/s,
                // so a ~100 s blend time-constant dominates it while staying
                // far from the acoustic dynamics. Using the full sponge
                // strength here (alpha ~0.9/step) shocked the dycore to NaN
                // within ~800 steps.
                Real alpha = amrex::min(Real(0.01) * xi * xi * dt, Real(0.1));
                amrex::ignore_unused(strength);
                Real rho_o = s_arr(i,j,k,Rho_comp);
                Real rho_n = (one - alpha) * rho_o + alpha * f_arr(i,j,k,Rho_comp);
                Real scale = rho_n / rho_o;
                s_arr(i,j,k,Rho_comp)       = rho_n;
                s_arr(i,j,k,RhoTheta_comp) *= scale;   // theta preserved
                if (blend_q) {
                    s_arr(i,j,k,RhoQ1_comp) *= scale;  // qv preserved
                }
            }
        });
    }
}

// Nest-band variant of the density blend (see nest-band comment in
// ERF::Advance). Free function for the same CUDA extended-lambda reason as
// hindcast_blend_band_density above.
static void nest_blend_band_density (MultiFab& S, const MultiFab& tgt,
                                     const iMultiFab& cf_mask, const int relax_val,
                                     const Real alpha, const bool has_moist)
{
    for (MFIter mfi(S); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const Array4<Real>& s_arr = S.array(mfi);
        const Array4<Real const>& t_arr = tgt.const_array(mfi);
        const Array4<int const>& m_arr = cf_mask.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (m_arr(i,j,k) == relax_val) {
                Real r_o = s_arr(i,j,k,Rho_comp);
                Real r_n = (one - alpha) * r_o + alpha * t_arr(i,j,k,Rho_comp);
                Real scale = r_n / r_o;
                s_arr(i,j,k,Rho_comp)       = r_n;
                s_arr(i,j,k,RhoTheta_comp) *= scale;
                if (has_moist) {
                    s_arr(i,j,k,RhoQ1_comp) *= scale;
                }
            }
        });
    }
}

static void bound_mynn_tke (MultiFab& S)
{
    const Real tke_max = Real(75.0);   // = qke_max/2, matching the internal clip
    for (MFIter mfi(S); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        const Array4<Real>& s_arr = S.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rho = s_arr(i,j,k,Rho_comp);
            s_arr(i,j,k,RhoKE_comp) = amrex::min(s_arr(i,j,k,RhoKE_comp), rho * tke_max);
        });
    }
}

/**
 * Function that advances the solution at one level for a single time step --
 * this does some preliminaries then calls erf_advance
 *
 * @param[in] lev level of refinement (coarsest level is 0)
 * @param[in] time start time for time advance
 * @param[in] dt_lev time step for this time advance
 */

void
ERF::Advance (int lev, Real time, Real dt_lev, int iteration, int /*ncycle*/)
{
    BL_PROFILE("ERF::Advance()");

    // We must swap the pointers so the previous step's "new" is now this step's "old"
    std::swap(vars_old[lev], vars_new[lev]);

    MultiFab& S_old = vars_old[lev][Vars::cons];
    MultiFab& S_new = vars_new[lev][Vars::cons];

    MultiFab& U_old = vars_old[lev][Vars::xvel];
    MultiFab& V_old = vars_old[lev][Vars::yvel];
    MultiFab& W_old = vars_old[lev][Vars::zvel];

    MultiFab& U_new = vars_new[lev][Vars::xvel];
    MultiFab& V_new = vars_new[lev][Vars::yvel];
    MultiFab& W_new = vars_new[lev][Vars::zvel];

    // We need to set these because otherwise in the first call to erf_advance we may
    //    read uninitialized data on ghost values in setting the bc's on the velocities
    U_new.setVal(Real(1.e34),U_new.nGrowVect());
    V_new.setVal(Real(1.e34),V_new.nGrowVect());
    W_new.setVal(Real(1.e34),W_new.nGrowVect());

    //
    // NOTE: the momenta here are not fillpatched (they are only used as scratch space)
    // If lev == 0 we have already FillPatched this in ERF::TimeStep
    //
    if (lev > 0) {
        // Set ghost cells to bogus values so they aren't uninitialized
        W_old.setBndry(Real(1.234e20));
        FillPatchFineLevel(lev, time, {&S_old, &U_old, &V_old, &W_old},
                           {&S_old, &rU_old[lev], &rV_old[lev], &rW_old[lev]},
                           base_state[lev], base_state[lev]);
    }

    // Momentum-consistent variant of the band-density blend (legacy call site
    // below, just before the dycore): running the blend BEFORE
    // VelocityToMomentum means the momenta are formed from the blended
    // density, so velocity is preserved exactly instead of being perturbed by
    // rho_old/rho_new in the band every step (UPSTREAM_ISSUES #12). The
    // FillBoundary syncs interior ghost rho so face averages at box seams see
    // blended values; physical-domain ghosts keep the pre-blend fill (the
    // outermost face factor is stale by O(alpha), and that cell is
    // relax-specified during the step anyway).
    static const bool blend_momentum = [] {
        bool b = false; ParmParse pp("erf"); pp.query("hindcast_blend_momentum", b); return b;
    }();
    static const bool blend_band_density_on = [] {
        bool b = false; ParmParse pp("erf"); pp.query("hindcast_blend_band_density", b); return b;
    }();
    if (blend_band_density_on && blend_momentum &&
        solverChoice.init_type == InitType::HindCast &&
        solverChoice.hindcast_lateral_forcing) {
        const bool l_has_moist = (solverChoice.moisture_type != MoistureType::None);
        hindcast_blend_band_density(S_old, forecast_state_interp[lev][Vars::cons],
                                    Geom(lev), dt_lev, solverChoice, l_has_moist);
        S_old.FillBoundary(Geom(lev).periodicity());
    }

    // Nest-band density blend (lev > 0): the c/f relaxation zone relaxes
    // momenta/theta/qv through the slow RHS but must NOT relax rho there
    // (mass-equation RHS relaxation destabilizes the acoustic substepping;
    // measured: SE-corner rho blowup in ~150 fine steps under the Jan-9
    // jet), yet with no rho constraint at all the band accumulates mass
    // (rho 1.18 -> 2.47 in ~450 fine steps, same run). Mirror the lev-0
    // solution: a gentle post-step blend of band rho toward the
    // time-interpolated parent, rescaling RhoTheta/RhoQ1 so theta and qv
    // are preserved, applied before VelocityToMomentum so momenta stay
    // consistent.
    static const Real l_cf_blend_alpha_gate = [] {
        Real a = Real(0.0); ParmParse pp("erf"); pp.query("cf_blend_alpha", a); return a;
    }();
    if (l_cf_blend_alpha_gate > Real(0.0) &&
        blend_band_density_on && lev > 0 && cf_width > 0 &&
        solverChoice.init_type == InitType::HindCast) {
        MultiFab tgt(S_old.boxArray(), S_old.DistributionMap(), S_old.nComp(), S_old.nGrowVect());
        MultiFab::Copy(tgt, S_old, 0, 0, S_old.nComp(), S_old.nGrowVect());
        PhysBCFunctNoOp void_bc;
        FPr_c[lev-1].FillRelax(tgt, time, void_bc, domain_bcs_type);
        auto* cf_mask = FPr_c[lev-1].GetMask();
        const int relax_val = FPr_c[lev-1].GetRelaxMaskVal();
        // Opt-in (erf.cf_blend_alpha > 0): the configuration verified on the
        // Sep-9 convective case uses rho RHS relaxation and NO nest blend;
        // this blend exists for the (still-open) strong-jet experiments
        // where rho relaxation is disabled via erf.cf_relax_rho=false.
        static const Real l_cf_blend_alpha = [] {
            Real a = Real(0.0); ParmParse pp("erf"); pp.query("cf_blend_alpha", a); return a;
        }();
        const Real alpha = std::min(l_cf_blend_alpha * dt_lev, Real(0.5));
        const bool l_has_moist = (solverChoice.moisture_type != MoistureType::None) &&
                                 (S_old.nComp() > RhoQ1_comp);
        if (alpha > Real(0.0)) {
            nest_blend_band_density(S_old, tgt, *cf_mask, relax_val, alpha, l_has_moist);
            S_old.FillBoundary(Geom(lev).periodicity());
        }
    }

    //
    // So we must convert the fillpatched to momenta, including the ghost values
    //
    const MultiFab* c_vfrac = nullptr;
    if (solverChoice.terrain_type == TerrainType::EB) {
        c_vfrac = &((get_eb(lev).get_const_factory())->getVolFrac());
    }

    VelocityToMomentum(U_old, rU_old[lev].nGrowVect(),
                       V_old, rV_old[lev].nGrowVect(),
                       W_old, rW_old[lev].nGrowVect(),
                       S_old, rU_old[lev], rV_old[lev], rW_old[lev],
                       Geom(lev).Domain(),
                       domain_bcs_type, c_vfrac);

    // Update the inflow perturbation update time and amplitude
    if (solverChoice.use_perturbation(lev))
    {
        turbPert.calc_tpi_update(lev, dt_lev, U_old, V_old, S_old);
    }

    // If PerturbationType::Direct or CPM is selected, directly add the computed perturbation
    // on the conserved field
    if (solverChoice.use_direct_perturbation(lev))
    {
        auto m_ixtype = S_old.boxArray().ixType(); // Conserved term
        for (MFIter mfi(S_old,TileNoZ()); mfi.isValid(); ++mfi) {
            Box bx  = mfi.tilebox();
            const Array4<Real> &cell_data  = S_old.array(mfi);
            const Array4<const Real> &pert_cell = turbPert.pb_cell[lev].array(mfi);
            turbPert.apply_tpi(lev, bx, RhoTheta_comp, m_ixtype, cell_data, pert_cell);
        }
    }

    // configure SurfaceLayer params if needed
    if (phys_bc_type[Orientation(Direction::z,Orientation::low)] == ERF_BC::surface_layer) {
        if (m_SurfaceLayer) {
            IntVect ng = Theta_prim[lev]->nGrowVect();
            MultiFab::Copy(  *Theta_prim[lev], S_old, RhoTheta_comp, 0, 1, ng);
            MultiFab::Divide(*Theta_prim[lev], S_old, Rho_comp     , 0, 1, ng);
            if (solverChoice.moisture_type != MoistureType::None) {
                ng = Qv_prim[lev]->nGrowVect();

                MultiFab::Copy(  *Qv_prim[lev], S_old, RhoQ1_comp, 0, 1, ng);
                MultiFab::Divide(*Qv_prim[lev], S_old, Rho_comp  , 0, 1, ng);

                if (solverChoice.moisture_indices.qr > -1) {
                    MultiFab::Copy(  *Qr_prim[lev], S_old, solverChoice.moisture_indices.qr, 0, 1, ng);
                    MultiFab::Divide(*Qr_prim[lev], S_old, Rho_comp  , 0, 1, ng);
                } else {
                    Qr_prim[lev]->setVal(0);
                }
            }
            // NOTE: std::swap above causes the field ptrs to be out of date.
            //       Reassign the field ptrs for MAC avg computation.
            m_SurfaceLayer->update_mac_ptrs(lev, vars_old, Theta_prim, Qv_prim, Qr_prim);
            m_SurfaceLayer->update_pblh(lev, vars_old, z_phys_cc[lev].get(),
                                        solverChoice.moisture_indices);

#ifdef ERF_USE_NETCDF
            Real elapsed_time_since_start_low = time + (start_time - start_low_time);
#else
            Real elapsed_time_since_start_low = time;
#endif
            m_SurfaceLayer->update_fluxes(lev, time, elapsed_time_since_start_low,
                                          S_old, z_phys_nd[lev], walldist[lev]);
        }
    }

#if defined(ERF_USE_WINDFARM)
    // **************************************************************************************
    // Update the windfarm sources
    // **************************************************************************************
    if (solverChoice.windfarm_type != WindFarmType::None) {
        advance_windfarm(Geom(lev), dt_lev, S_old,
                         U_old, V_old, W_old, vars_windfarm[lev],
                         Nturb[lev], SMark[lev], time);
    }

#endif

    // **************************************************************************************
    // Update the radiation sources with the "old" state
    // **************************************************************************************
    advance_radiation(lev, S_old, dt_lev);

#ifdef ERF_USE_SHOC
    // **************************************************************************************
    // Update the "old" state using SHOC
    // **************************************************************************************
    if (solverChoice.use_shoc) {
        // Get SFC fluxes from SurfaceLayer
        if (m_SurfaceLayer) {
            Vector<const MultiFab*> mfs = {&S_old, &U_old, &V_old, &W_old};
            m_SurfaceLayer->impose_SurfaceLayer_bcs(lev, mfs, Tau[lev],
                                                    SFS_hfx1_lev[lev].get() , SFS_hfx2_lev[lev].get() , SFS_hfx3_lev[lev].get(),
                                                    SFS_q1fx1_lev[lev].get(), SFS_q1fx2_lev[lev].get(), SFS_q1fx3_lev[lev].get(),
                                                    z_phys_nd[lev].get());
        }

        // Get Shoc tendencies and update the state
        Real* w_sub = (solverChoice.custom_w_subsidence) ? d_w_subsid[lev].data() : nullptr;
        compute_shoc_tendencies(lev, &S_old, &U_old, &V_old, &W_old, w_sub,
                                Tau[lev][TauType::tau13].get(), Tau[lev][TauType::tau23].get(),
                                SFS_hfx3_lev[lev].get()       , SFS_q1fx3_lev[lev].get()      ,
                                eddyDiffs_lev[lev].get()      , z_phys_nd[lev].get()          ,
                                dt_lev);
    }
#endif

    const BoxArray&            ba = S_old.boxArray();
    const DistributionMapping& dm = S_old.DistributionMap();

    int nvars = S_old.nComp();

    // Source array for conserved cell-centered quantities -- this will be filled
    //     in the call to make_sources in ERF_TI_slow_rhs_pre.H
    MultiFab cc_source(ba,dm,nvars,1); cc_source.setVal(0);

    // Source arrays for momenta -- these will be filled
    //     in the call to make_mom_sources in ERF_TI_slow_rhs_pre.H
    BoxArray ba_x(ba); ba_x.surroundingNodes(0);
    MultiFab xmom_source(ba_x,dm,1,1); xmom_source.setVal(0);

    BoxArray ba_y(ba); ba_y.surroundingNodes(1);
    MultiFab ymom_source(ba_y,dm,1,1); ymom_source.setVal(0);

    BoxArray ba_z(ba); ba_z.surroundingNodes(2);
    MultiFab zmom_source(ba_z,dm,1,1); zmom_source.setVal(0);
    MultiFab    buoyancy(ba_z,dm,1,1); buoyancy.setVal(0);

    amrex::Vector<MultiFab> state_old;
    amrex::Vector<MultiFab> state_new;

    // **************************************************************************************
    // Here we define state_old and state_new which are to be advanced
    // **************************************************************************************
    // Initial solution
    // Note that "old" and "new" here are relative to each RK stage.
    state_old.push_back(MultiFab(S_old      , amrex::make_alias, 0, nvars)); // cons
    state_old.push_back(MultiFab(rU_old[lev], amrex::make_alias, 0,     1)); // xmom
    state_old.push_back(MultiFab(rV_old[lev], amrex::make_alias, 0,     1)); // ymom
    state_old.push_back(MultiFab(rW_old[lev], amrex::make_alias, 0,     1)); // zmom

    // Final solution
    // state_new at the end of the last RK stage holds the t^{n+1} data
    state_new.push_back(MultiFab(S_new      , amrex::make_alias, 0, nvars)); // cons
    state_new.push_back(MultiFab(rU_new[lev], amrex::make_alias, 0,     1)); // xmom
    state_new.push_back(MultiFab(rV_new[lev], amrex::make_alias, 0,     1)); // ymom
    state_new.push_back(MultiFab(rW_new[lev], amrex::make_alias, 0,     1)); // zmom

    // **************************************************************************************
    // Bound the prognostic MYNN TKE like WRF does. The MYNN-EDMF column solver
    // clips its internal qke to <= 150 m2/s2 (2*TKE), and in WRF that clipped
    // qke IS the prognostic state; in this port the dycore-advected RhoKE never
    // sees the clip, and upper-level shear production (2.5-km cells at 8-15 km)
    // was measured running TKE past 200 m2/s2 to NaN at ~4.6 h. Apply the same
    // bound to the state.
    // **************************************************************************************
    static const bool bound_pbl_tke = [] {
        bool b = true; ParmParse pp("erf"); pp.query("bound_pbl_tke", b); return b;
    }();
    if (bound_pbl_tke &&
        (solverChoice.turbChoice[lev].pbl_type == PBLType::MYNNEDMF ||
         solverChoice.turbChoice[lev].pbl_type == PBLType::MYNN25)) {
        bound_mynn_tke(S_old);
    }

    // Direct band-density blending: three rates were tried (0.9/step, 0.04/step,
    // and an RHS mass source) and ALL destabilize the acoustic dycore faster
    // than the slow band mass-accumulation they target. Default OFF; kept for
    // future specified-zone boundary work.
    // (Skipped when erf.hindcast_blend_momentum is set: the momentum-consistent
    // variant already ran before VelocityToMomentum above.)
    static const bool blend_band_density = [] {
        bool b = false; ParmParse pp("erf"); pp.query("hindcast_blend_band_density", b); return b;
    }();
    if (blend_band_density && !blend_momentum &&
        solverChoice.init_type == InitType::HindCast &&
        solverChoice.hindcast_lateral_forcing) {
        const bool l_has_moist = (solverChoice.moisture_type != MoistureType::None);
        hindcast_blend_band_density(S_old, forecast_state_interp[lev][Vars::cons],
                                    Geom(lev), dt_lev, solverChoice, l_has_moist);
    }

    // **************************************************************************************
    // Tests on the reasonableness of the solution before the dycore
    // **************************************************************************************
    // Test for NaNs after dycore
    if (check_for_nans > 1 && nan_check_step(istep[0])) {
        if (verbose > 1) {
            amrex::Print() << "Testing old state and vels for NaNs before dycore" << std::endl;
        }
        check_state_for_nans(S_old);
        check_vels_for_nans(rU_old[lev],rV_old[lev],rW_old[lev]);
    }

    // We only test on low temp if we have a moisture model because we are protecting against
    //    the test on low temp inside the moisture models
    if (solverChoice.moisture_type != MoistureType::None) {
        if (verbose > 1) {
            amrex::Print() << "Testing on low temperature before dycore" << std::endl;
        }
        check_for_low_temp(S_old);
    } else {
        if (verbose > 1) {
            amrex::Print() << "Testing on negative temperature before dycore" << std::endl;
        }
        check_for_negative_theta(S_old);
    }

    debug_cell_print(S_old, "pre_dycore ");

    // **************************************************************************************
    // Update the dycore
    // **************************************************************************************
    advance_dycore(lev, state_old, state_new,
                   U_old, V_old, W_old,
                   U_new, V_new, W_new,
                   cc_source, xmom_source, ymom_source, zmom_source, buoyancy,
                   Geom(lev), dt_lev, time);

    // **************************************************************************************
    // Tests on the reasonableness of the solution after the dycore
    // **************************************************************************************
    // Test for NaNs after dycore
    if (check_for_nans > 0 && nan_check_step(istep[0])) {
        if (verbose > 1) {
            amrex::Print() << "Testing new state and vels for NaNs after dycore" << std::endl;
        }
        check_state_for_nans(S_new);
        check_vels_for_nans(rU_new[lev],rV_new[lev],rW_new[lev]);
    }

    // We only test on low temp if we have a moisture model because we are protecting against
    //    the test on low temp inside the moisture models
    if (solverChoice.moisture_type != MoistureType::None) {
        if (verbose > 1) {
            amrex::Print() << "Testing on low temperature after dycore" << std::endl;
        }
        check_for_low_temp(S_new);
    } else {
        // Otherwise we will test on negative (rhotheta) coming out of the dycore
        if (verbose > 1) {
            amrex::Print() << "Testing on negative temperature after dycore" << std::endl;
        }
        check_for_negative_theta(S_new);
    }

    debug_cell_print(S_new, "post_dycore");

    // **************************************************************************************
    // Update the microphysics (moisture)
    // **************************************************************************************
    if (!solverChoice.moisture_tight_coupling)
    {
        advance_microphysics(lev, S_new, dt_lev, iteration, time);

        // Test for NaNs after microphysics
        if (check_for_nans > 0 && nan_check_step(istep[0])) {
            amrex::Print() << "Testing new state for NaNs after advance_microphysics" << std::endl;
            check_state_for_nans(S_new);
        }
    }

    debug_cell_print(S_new, "post_micro ");

    // **************************************************************************************
    // Update the land surface model
    // **************************************************************************************
    Real time_at_end_of_step = time+dt_lev;
    advance_lsm(lev, S_new, U_new, V_new, time_at_end_of_step, dt_lev);

#ifdef ERF_USE_PARTICLES
    // **************************************************************************************
    // Update the particle positions
    // **************************************************************************************
   evolveTracers(lev, dt_lev, vars_new, z_phys_nd);
#endif

    // ***********************************************************************************************
    // Impose domain boundary conditions here so that in FillPatching the fine data we won't
    // need to re-fill these
    // ***********************************************************************************************
    if (lev < finest_level) {
         IntVect ngvect_vels = vars_new[lev][Vars::xvel].nGrowVect();
         (*physbcs_cons[lev])(vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              0,vars_new[lev][Vars::cons].nComp(),
                              vars_new[lev][Vars::cons].nGrowVect(),time,BCVars::cons_bc,true);
            (*physbcs_u[lev])(vars_new[lev][Vars::xvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::xvel_bc,true);
            (*physbcs_v[lev])(vars_new[lev][Vars::yvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::yvel_bc,true);
            (*physbcs_w[lev])(vars_new[lev][Vars::zvel], vars_new[lev][Vars::xvel], vars_new[lev][Vars::yvel],
                              ngvect_vels,time,BCVars::zvel_bc,true);
    }

    // **************************************************************************************
    // Register old and new coarse data if we are at a level less than the finest level
    // **************************************************************************************
    if (lev < finest_level) {
        if (cf_width > 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::cons].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::cons].FillBoundary(geom[lev].periodicity());
            FPr_c[lev].RegisterCoarseData({&state_old[IntVars::cons], &state_new[IntVars::cons]},
                                          {time, time+dt_lev});
        }

        if (cf_width >= 0) {
            // We must fill the ghost cells of these so that the parallel copy works correctly
            state_old[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::xmom].FillBoundary(geom[lev].periodicity());
            FPr_u[lev].RegisterCoarseData({&state_old[IntVars::xmom], &state_new[IntVars::xmom]},
                                          {time, time+dt_lev});

            state_old[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::ymom].FillBoundary(geom[lev].periodicity());
            FPr_v[lev].RegisterCoarseData({&state_old[IntVars::ymom], &state_new[IntVars::ymom]},
                                          {time, time+dt_lev});

            state_old[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            state_new[IntVars::zmom].FillBoundary(geom[lev].periodicity());
            FPr_w[lev].RegisterCoarseData({&state_old[IntVars::zmom], &state_new[IntVars::zmom]},
                                          {time, time+dt_lev});
        }

            //
            // Now create a MultiFab that holds (S_new - S_old) / dt from the coarse level interpolated
            //     on to the coarse/fine boundary at the fine resolution
            //
            Interpolater* mapper_f = &face_cons_linear_interp;

            // PhysBCFunctNoOp null_bc;
            // MultiFab tempx(vars_new[lev+1][Vars::xvel].boxArray(),vars_new[lev+1][Vars::xvel].DistributionMap(),1,0);
            // tempx.setVal(0);
            // xmom_crse_rhs[lev+1].setVal(0);
            // FPr_u[lev].FillSet(tempx               , time       , null_bc, domain_bcs_type);
            // FPr_u[lev].FillSet(xmom_crse_rhs[lev+1], time+dt_lev, null_bc, domain_bcs_type);
            // MultiFab::Subtract(xmom_crse_rhs[lev+1],tempx,0,0,1,IntVect{0});
            // xmom_crse_rhs[lev+1].mult(one/dt_lev,0,1,0);

            // MultiFab tempy(vars_new[lev+1][Vars::yvel].boxArray(),vars_new[lev+1][Vars::yvel].DistributionMap(),1,0);
            // tempy.setVal(0);
            // ymom_crse_rhs[lev+1].setVal(0);
            // FPr_v[lev].FillSet(tempy               , time       , null_bc, domain_bcs_type);
            // FPr_v[lev].FillSet(ymom_crse_rhs[lev+1], time+dt_lev, null_bc, domain_bcs_type);
            // MultiFab::Subtract(ymom_crse_rhs[lev+1],tempy,0,0,1,IntVect{0});
            // ymom_crse_rhs[lev+1].mult(one/dt_lev,0,1,0);

            MultiFab temp_state(zmom_crse_rhs[lev+1].boxArray(),zmom_crse_rhs[lev+1].DistributionMap(),1,0);
            InterpFromCoarseLevel(temp_state,            IntVect{0}, IntVect{0}, state_old[IntVars::zmom], 0, 0, 1,
                                  geom[lev], geom[lev+1], refRatio(lev), mapper_f, domain_bcs_type, BCVars::zvel_bc);
            InterpFromCoarseLevel(zmom_crse_rhs[lev+1],  IntVect{0}, IntVect{0}, state_new[IntVars::zmom], 0, 0, 1,
                                  geom[lev], geom[lev+1], refRatio(lev), mapper_f, domain_bcs_type, BCVars::zvel_bc);
            MultiFab::Subtract(zmom_crse_rhs[lev+1],temp_state,0,0,1,IntVect{0});
            zmom_crse_rhs[lev+1].mult(one/dt_lev,0,1,0);
    }

    // ***********************************************************************************************
    // Update the time averaged velocities if they are requested
    // ***********************************************************************************************
    if (solverChoice.time_avg_vel) {
        Time_Avg_Vel_atCC(dt[lev], t_avg_cnt[lev], vel_t_avg[lev].get(), U_new, V_new, W_new);
    }
}
