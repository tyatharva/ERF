#include <ERF.H>

using namespace amrex;

void ERF::advance_radiation (int lev,
                             MultiFab& cons,
                             const Real& dt_advance)
{
    if (solverChoice.rad_type != RadiationType::None) {
#ifdef ERF_USE_NETCDF
        MultiFab *lat_ptr = lat_m[lev].get();
        MultiFab *lon_ptr = lon_m[lev].get();
#else
        MultiFab *lat_ptr = nullptr;
        MultiFab *lon_ptr = nullptr;
#endif
        // T surf from SurfaceLayer if we have it
        MultiFab* t_surf = (m_SurfaceLayer) ? m_SurfaceLayer->get_t_surf(lev) : nullptr;

        // RRTMGP inputs names and pointers
        Vector<std::string> lsm_input_names = rad[lev]->get_lsm_input_varnames();
        Vector<MultiFab*> lsm_input_ptrs(lsm_input_names.size(),nullptr);
        for (int i(0); i<lsm_input_ptrs.size(); ++i) {
            int varIdx = lsm.Get_DataIdx(lev,lsm_input_names[i]);
            if (varIdx >= 0) { lsm_input_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
        }

        // RRTMGP output names and pointers
        Vector<std::string> lsm_output_names = rad[lev]->get_lsm_output_varnames();
        Vector<MultiFab*> lsm_output_ptrs(lsm_output_names.size(),nullptr);
        for (int i(0); i<lsm_output_ptrs.size(); ++i) {
            int varIdx = lsm.Get_DataIdx(lev,lsm_output_names[i]);
            if (varIdx >= 0) { lsm_output_ptrs[i] = lsm.Get_Data_Ptr(lev,varIdx); }
        }

        // Anelastic: radiation thermodynamics must use the hydrostatic
        // reference pressure p0(z), not EOS of the frozen-density state
        // (see Radiation::set_base_state).
        if (solverChoice.anelastic[lev]) {
            rad[lev]->set_base_state(&base_state[lev]);
        }

        // Per-column surface albedo from the hindcast surface frames
        // (ERA5 fal), when the frames carry it.
        if (!alb_lev[lev].empty() && alb_lev[lev][0]) {
            rad[lev]->set_albedo(alb_lev[lev][0].get());
        }

        // The rad_freq_in_time trigger runs on MODEL time, not on the absolute
        // instant below: t_old + start_time is a float32 sum of a ~1.673e9
        // epoch value (ULP 128 s), which quantised the cadence to multiples of
        // 128 s and turned a 180 s request into a measured 255.4 s
        // (UPSTREAM_ISSUES #26, the same float32-epoch class as #28).
        rad[lev]->set_model_time(t_old[lev]);

        // Enter radiation class driver
        amrex::Real time_for_rad = t_old[lev] + start_time;
        rad[lev]->Run(lev, istep[lev], time_for_rad, dt_advance,
                      cons.boxArray(), geom[lev], &(cons),
                      lmask_lev[lev][0].get(), t_surf,
                      lsm_input_ptrs, lsm_output_ptrs,
                      qheating_rates[lev].get(), rad_fluxes[lev].get(),
                      z_phys_nd[lev].get()     , lat_ptr, lon_ptr);
    }
}
