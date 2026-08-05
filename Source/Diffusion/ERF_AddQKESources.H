    // Using PBL
    if (l_use_mynn && start_comp <= RhoKE_comp && end_comp >=RhoKE_comp) {
        int qty_index = RhoKE_comp;
        auto pbl_mynn_B1_l = turbChoice.pbl_mynn.B1;

        // Surface TKE production. See ERF_TurbStruct.H for why this is needed at
        // all; the form below is transcribed from the WRF-MYNN port that sits
        // unused behind `#if 0` in ERF_ComputeDiffusivityMYNNEDMF.cpp:1097,1130:
        //     vkz  = karman * 0.5 * dz[kts]
        //     pdk1 = 2 * ust^3 * pmz / vkz
        // pmz is phi_m - zeta evaluated at z1, and defaults to 1 in that port
        // (line 4079); we keep that default rather than invent a stability
        // correction. WRF's pdk1 is a source for QKE = 2*TKE, and this routine
        // supplies a source for TKE (see the note in ComputeQKESourceTerms:
        // "we transport TKE = myhalf*QKE ... so the RHS terms do not have a
        // factor of two"), so the factor of 2 is dropped here and the term is
        //     rho * ust^3 / (karman * z1),    z1 = 0.5*dz
        // Note z1 uses the NOMINAL dz, not the terrain-following cell thickness
        // (19.7-30.7 m against a nominal 25 m on Domain A). That is the same
        // approximation the shear-production term above already makes -- this
        // call passes no met_h_zeta, so it defaults to 1 -- and being more
        // precise in one term than its neighbours would not be an improvement.
        const bool l_surf_tke_src = turbChoice.pbl_surface_tke_source;
        const Real karman_l  = KAPPA;
        const Real z1_l      = Real(0.5) / cellSizeInv[2];
        const int  klo_l     = domain.smallEnd(2);

        ParallelFor(bx,[=,moist_indices=solverChoice.moisture_indices] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
            bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].hi(2) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
            bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir) );

            // This computes shear production, buoyancy production, and dissipation terms only.
            cell_rhs(i, j, k, qty_index) += ComputeQKESourceTerms(i,j,k,u,v,cell_data,cell_prim,
                                                                  mu_turb,cellSizeInv,domain,
                                                                  pbl_mynn_B1_l,tm_arr(i,j,0),
                                                                  moist_indices,
                                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                                                  u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                                                  v_ext_dir_on_zlo, v_ext_dir_on_zhi);

            if (l_surf_tke_src && k == klo_l && u_star_arr) {
                const Real ust = u_star_arr(i,j,0);
                cell_rhs(i, j, k, qty_index) += cell_data(i,j,k,Rho_comp)
                                              * ust*ust*ust / (karman_l * z1_l);
            }
        });
    }
