#ifndef ERF_WEATHERDATAINTERPOLATION_H_
#define ERF_WEATHERDATAINTERPOLATION_H_
/**
 * Trilinear interpolation of weather forecast data onto the simulation mesh
 * The coarse weather forecast data is interpolated in time first to get the forecast
 * at the current time, and then spatially interpolated onto the simulation mesh
 */

#include <filesystem>
#include <stdexcept>
#include "ERF.H"
#include "ERF_ReadCustomBinaryIC.H"
#include "ERF_Interpolation_Bilinear.H"

using namespace amrex;
namespace fs = std::filesystem;

enum class MultiFabType { CC, NC };

void PlotMultiFab(const MultiFab& mf,
                  const Geometry& geom_mf,
                  const std::string plotfilename,
                  MultiFabType mftype)
{

    Vector<std::string> varnames = {
    "rho", "uvel", "vvel", "wvel", "theta", "qv", "qc", "qr", "latitude", "longitude"
    }; // Customize variable names

    const Real time = zero;


    // Assume weather_mf is nodal in all directions
    if(mftype == MultiFabType::NC) {
        BoxArray cba = mf.boxArray();
        cba = amrex::convert(mf.boxArray(), IntVect::TheCellVector());

        MultiFab cc_mf(cba, mf.DistributionMap(),
               mf.nComp(), 0);

        amrex::average_node_to_cellcenter(cc_mf, 0, mf, 0, mf.nComp());

        WriteSingleLevelPlotfile(
            plotfilename,
            cc_mf,
            varnames,
            geom_mf,
            time,
            0 // level
        );
    } else {
        WriteSingleLevelPlotfile(
            plotfilename,
            mf,
            varnames,
            geom_mf,
            time,
            0 // level
        );
    }
}

void
ERF::FillForecastStateMultiFabs(const int lev,
                                const std::string& filename,
                                const std::unique_ptr<MultiFab>& a_z_phys_nd,
                                Vector<Vector<MultiFab>>& forecast_state)
{

    Vector<Real> latvec_h, lonvec_h, xvec_h, yvec_h, zvec_h;
    Vector<Real> rho_h, uvel_h, vvel_h, wvel_h, theta_h, qv_h, qc_h, qr_h;

    ReadCustomBinaryIC(filename, latvec_h, lonvec_h,
                       xvec_h, yvec_h, zvec_h, rho_h,
                       uvel_h, vvel_h, wvel_h,
                       theta_h, qv_h, qc_h, qr_h);

    Real zmax = *std::max_element(zvec_h.begin(), zvec_h.end());

    const auto prob_lo_erf  = geom[lev].ProbLoArray();
    const auto prob_hi_erf  = geom[lev].ProbHiArray();
    const auto dx_erf       = geom[lev].CellSizeArray();

    if (prob_hi_erf[2] >= zmax) {
        Abort("ERROR: the maximum z of the domain (" + std::to_string(prob_hi_erf[2]) +
        ") should be less than the maximum z in the forecast data (" + std::to_string(zmax) +
        "). Change geometry.prob_hi[2] in the inputs to be less than " + std::to_string(zmax) + "."
        );
    }

    if(prob_lo_erf[0] < xvec_h.front() + 4*dx_erf[0]){
        amrex::Abort("The xlo value of the domain has to be greater than " + std::to_string(xvec_h.front() + 4*dx_erf[0]));
    }
    if(prob_hi_erf[0] > xvec_h.back() - 4*dx_erf[0]){
        amrex::Abort("The xhi value of the domain has to be less than " + std::to_string(xvec_h.back() - 4*dx_erf[0]));
    }
    if(prob_lo_erf[1] < yvec_h.front() + 4*dx_erf[1]){
        amrex::Abort("The ylo value of the domain has to be greater than " + std::to_string(yvec_h.front() + 4*dx_erf[1]));
    }
    if(prob_hi_erf[1] > yvec_h.back() - 4*dx_erf[1]){
        amrex::Abort("The yhi value of the domain has to be less than " + std::to_string(yvec_h.back() - 4*dx_erf[1]));
    }


    int nx = static_cast<int>(xvec_h.size());
    int ny = static_cast<int>(yvec_h.size());
    int nz = static_cast<int>(zvec_h.size());

    amrex::Real dxvec = (xvec_h[nx-1]-xvec_h[0])/(nx-1);
    amrex::Real dyvec = (yvec_h[ny-1]-yvec_h[0])/(ny-1);

    amrex::Gpu::DeviceVector<Real> latvec_d(nx*ny), lonvec_d(nx*ny), zvec_d(nz);
    amrex::Gpu::DeviceVector<Real> xvec_d(nx*ny*nz), yvec_d(nx*ny*nz);
    amrex::Gpu::DeviceVector<Real> rho_d(nx*ny*nz), uvel_d(nx*ny*nz), vvel_d(nx*ny*nz), wvel_d(nx*ny*nz),
                                   theta_d(nx*ny*nz), qv_d(nx*ny*nz), qc_d(nx*ny*nz), qr_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, latvec_h.begin(), latvec_h.end(), latvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, lonvec_h.begin(), lonvec_h.end(), lonvec_d.begin());

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xvec_h.begin(), xvec_h.end(), xvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yvec_h.begin(), yvec_h.end(), yvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, rho_h.begin(), rho_h.end(), rho_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, theta_h.begin(), theta_h.end(), theta_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, uvel_h.begin(), uvel_h.end(), uvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, vvel_h.begin(), vvel_h.end(), vvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, wvel_h.begin(), wvel_h.end(), wvel_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qv_h.begin(), qv_h.end(), qv_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qc_h.begin(), qc_h.end(), qc_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, qr_h.begin(), qr_h.end(), qr_d.begin());

    amrex::Gpu::streamSynchronize();

    Real* latvec_d_ptr = latvec_d.data();
    Real* lonvec_d_ptr = lonvec_d.data();
    Real* xvec_d_ptr = xvec_d.data();
    Real* yvec_d_ptr = yvec_d.data();
    Real* zvec_d_ptr = zvec_d.data();
    Real* rho_d_ptr   = rho_d.data();
    Real* uvel_d_ptr  = uvel_d.data();
    Real* vvel_d_ptr  = vvel_d.data();
    Real* wvel_d_ptr  = wvel_d.data();
    Real* theta_d_ptr = theta_d.data();
    Real* qv_d_ptr = qv_d.data();
    Real* qc_d_ptr = qc_d.data();
    Real* qr_d_ptr = qr_d.data();

    MultiFab& erf_mf_cons   = forecast_state[lev][Vars::cons];
    MultiFab& erf_mf_xvel   = forecast_state[lev][Vars::xvel];
    MultiFab& erf_mf_yvel   = forecast_state[lev][Vars::yvel];
    MultiFab& erf_mf_zvel   = forecast_state[lev][Vars::zvel];
    MultiFab& erf_mf_latlon = forecast_state[lev][4];

    erf_mf_cons.setVal(0.0);
    erf_mf_xvel.setVal(0.0);
    erf_mf_yvel.setVal(0.0);
    erf_mf_zvel.setVal(0.0);
    erf_mf_latlon.setVal(0.0);

    // Interpolate the data on to the ERF mesh

     for (MFIter mfi(erf_mf_cons); mfi.isValid(); ++mfi) {
        const auto z_arr    = (a_z_phys_nd) ? a_z_phys_nd->const_array(mfi) :
                                            Array4<const Real> {};
        const Array4<Real> &fine_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &fine_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &fine_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &fine_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &fine_latlon_arr = erf_mf_latlon.array(mfi);


        const Box& gbx = mfi.growntilebox(); // tilebox + ghost cells
        const int ncomp_cons = erf_mf_cons.nComp();

        const Box &gtbx = mfi.tilebox(IntVect(1,0,0));
        const Box &gtby = mfi.tilebox(IntVect(0,1,0));
        const Box &gtbz = mfi.tilebox(IntVect(0,0,1));
        const auto prob_lo  = geom[lev].ProbLoArray();
        const auto dx       = geom[lev].CellSizeArray();
       //const Box &gtbz = mfi.tilebox(IntVect(0,0,1));

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Geometry (note we must include these here to get the data on device)
            const Real x        = prob_lo[0] + (i + myhalf) * dx[0];
            const Real y        = prob_lo[1] + (j + myhalf) * dx[1];
            //const Real z        = prob_lo[2] + (k + myhalf) * dx[2];
            const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/two;

            // First interpolate where the weather data is available from
            Real tmp_rho, tmp_theta, tmp_qv, tmp_qc, tmp_qr, tmp_lat, tmp_lon;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   rho_d_ptr, tmp_rho);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   theta_d_ptr, tmp_theta);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   qv_d_ptr, tmp_qv);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   qc_d_ptr, tmp_qc);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   qr_d_ptr, tmp_qr);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, 1,
                                   x, y, zero,
                                   latvec_d_ptr, tmp_lat);

            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, 1,
                                   x, y, zero,
                                   lonvec_d_ptr, tmp_lon);

            fine_cons_arr(i,j,k,Rho_comp) = tmp_rho;
            // Store PLAIN theta and qv (not rho-weighted), mirroring the
            // velocity convention: the boundary sponge forms rho_f*theta_f /
            // rho_f*qv_f at application time. Slots beyond Rho_comp were
            // previously discarded (momenta-only coupling gap).
            fine_cons_arr(i,j,k,RhoTheta_comp) = tmp_theta;
            if (ncomp_cons > RhoQ1_comp) {
                fine_cons_arr(i,j,k,RhoQ1_comp) = tmp_qv;
            }
            fine_latlon_arr(i,j,k,0) = tmp_lat;
            fine_latlon_arr(i,j,k,1) = tmp_lon;
        });

        ParallelFor(gtbx, gtby, gtbz,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + i       * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+myhalf) * dx_erf[1];
            //Real z = prob_lo_erf[2] + (k+myhalf) * dx_erf[2];
            const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/two;

            Real tmp_uvel;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   uvel_d_ptr, tmp_uvel);

            fine_xvel_arr(i, j, k, 0) = tmp_uvel;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+myhalf) * dx_erf[0];
            Real y = prob_lo_erf[1] + j       * dx_erf[1];
            //Real z = prob_lo_erf[2] + (k+myhalf) * dx_erf[2];
            const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/two;

            Real tmp_vvel;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   vvel_d_ptr, tmp_vvel);

            fine_yvel_arr(i, j, k, 0) = tmp_vvel;
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
             // Physical location of the fine node
            Real x = prob_lo_erf[0] + (i+myhalf) * dx_erf[0];
            Real y = prob_lo_erf[1] + (j+myhalf) * dx_erf[1];
            Real z = prob_lo_erf[2] + k       * dx_erf[2];
            //const Real z = (z_arr(i,j,k) + z_arr(i,j,k+1))/two;

            Real tmp_wvel;
            bilinear_interpolation(xvec_d_ptr, yvec_d_ptr, zvec_d_ptr,
                                   dxvec, dyvec,
                                   nx, ny, nz,
                                   x, y, z,
                                   wvel_d_ptr, tmp_wvel);

            fine_zvel_arr(i, j, k, 0) = tmp_wvel;
        });
    }

    /*Vector<std::string> varnames = {
    "rho", "uvel", "vvel", "wvel", "theta", "qv", "qc", "qr"
    }; // Customize variable names

     Vector<std::string> varnames_cons = {
    "rho", "rhotheta", "ke", "sc", "rhoqv", "rhoqc", "rhoqr"
    }; // Customize variable names

    Vector<std::string> varnames_plot_mf = {
    "rho", "rhotheta", "rhoqv", "rhoqc", "rhoqr", "xvel", "yvel", "zvel", "latitude", "longitude"
    }; // Customize variable names

    const Real time = zero;

    std::string pltname = "plt_interp";

    MultiFab plot_mf(erf_mf_cons.boxArray(), erf_mf_cons.DistributionMap(),
                     10, 0);

    plot_mf.setVal(0.0);

    for (MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
        const Array4<Real> &plot_mf_arr = plot_mf.array(mfi);
        const Array4<Real> &erf_mf_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &erf_mf_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &erf_mf_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &erf_mf_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &erf_mf_latlon_arr = erf_mf_latlon.array(mfi);

        const Box& bx = mfi.validbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plot_mf_arr(i,j,k,0) = erf_mf_cons_arr(i,j,k,Rho_comp);
            plot_mf_arr(i,j,k,1) = erf_mf_cons_arr(i,j,k,RhoTheta_comp);
            plot_mf_arr(i,j,k,2) = erf_mf_cons_arr(i,j,k,RhoQ1_comp);
            plot_mf_arr(i,j,k,3) = erf_mf_cons_arr(i,j,k,RhoQ2_comp);
            plot_mf_arr(i,j,k,4) = erf_mf_cons_arr(i,j,k,RhoQ3_comp);

            plot_mf_arr(i,j,k,5) = (erf_mf_xvel_arr(i,j,k,0) + erf_mf_xvel_arr(i+1,j,k,0))/two;
            plot_mf_arr(i,j,k,6) = (erf_mf_yvel_arr(i,j,k,0) + erf_mf_yvel_arr(i,j+1,k,0))/two;
            plot_mf_arr(i,j,k,7) = (erf_mf_zvel_arr(i,j,k,0) + erf_mf_zvel_arr(i,j,k+1,0))/two;

            plot_mf_arr(i,j,k,8) = erf_mf_latlon_arr(i,j,k,0);
            plot_mf_arr(i,j,k,9) = erf_mf_latlon_arr(i,j,k,1);
        });
    }


    WriteSingleLevelPlotfile(
            pltname,
            plot_mf,
            varnames_plot_mf,
            geom[0],
            time,
            0 // level
        );*/
}

#ifdef ERF_USE_NETCDF
/**
 * Copy one component of a (possibly staggered) MultiFab into a global
 * FArrayBox covering `strip` that is identical on all ranks: ParallelCopy
 * gathers the distributed data onto a single-box MultiFab on rank 0, which
 * is then broadcast. (Same all-ranks convention as the metgrid/wrfbdy
 * boundary planes; ParallelCopy also resolves shared staggered faces
 * without double counting.)
 */
static void
strip_to_global_fab (const MultiFab& src, const int scomp,
                     const Box& strip, FArrayBox& dest)
{
    BoxArray sba(strip);
    DistributionMapping sdm(Vector<int>({0}));
    MultiFab smf(sba, sdm, 1, 0);
    smf.ParallelCopy(src, scomp, 0, 1);

    const Long npts = strip.numPts();
    Vector<Real> buf(npts, Real(0.0));
    for (MFIter mfi(smf); mfi.isValid(); ++mfi) {   // non-empty on rank 0 only
        Gpu::copy(Gpu::deviceToHost,
                  smf[mfi].dataPtr(), smf[mfi].dataPtr() + npts, buf.data());
    }
    ParallelDescriptor::Bcast(buf.data(), npts, 0);
    Gpu::copy(Gpu::hostToDevice, buf.data(), buf.data() + npts, dest.dataPtr());
}

/**
 * Fill the real-BC lateral boundary planes from the hindcast frames.
 *
 * The specified+relaxation-zone machinery (fill_from_realbdy each FillPatch;
 * realbdy_compute_interior_ghost_rhs in the slow RHS) consumes
 * bdy_data_{xlo,xhi,ylo,yhi}[time][RealBdyVars::{U,V,T,QV}] as global
 * FArrayBox strips of width real_width holding PLAIN u, v, theta, qv --
 * exactly what FillForecastStateMultiFabs interpolates from the ERA5/GFS
 * .bin frames. This routine sizes the planes like init_from_metgrid does
 * and fills every frame, so init_type = HindCast gets the WRF-style inflow
 * treatment (boundary values SET, not just interior-relaxed).
 */
void
ERF::fill_bdy_data_from_hindcast ()
{
    const int lev = 0;

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(real_width > 0,
        "erf.use_real_bcs with HindCast requires erf.real_width > 0");

    // Enumerate the frames the same way WeatherDataInterpolation does
    std::string folder = solverChoice.hindcast_boundary_data_dir;
    std::vector<std::string> bin_files;
    for (const auto& entry : fs::directory_iterator(folder)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        if (fname.size() >= 4 && fname.substr(fname.size() - 4) == ".bin") {
            bin_files.push_back(entry.path().string());
        }
    }
    std::sort(bin_files.begin(), bin_files.end());
    const int ntimes = static_cast<int>(bin_files.size());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ntimes >= 2,
        "Need at least two hindcast frames for boundary data");

    // The first .bin frame corresponds to start_datetime by construction
    // (WeatherDataInterpolation indexes frames by elapsed time/interval).
    bdy_time_interval = solverChoice.hindcast_data_interval_in_hrs * Real(3600.0);
    start_bdy_time    = start_time;
    final_bdy_time    = start_time + (ntimes-1) * bdy_time_interval;

    const bool l_use_moisture = (solverChoice.moisture_type != MoistureType::None);
    const int BdyEnd0 = l_use_moisture ? MetGridBdyVars::NumTypes
                                       : MetGridBdyVars::NumTypes-1;

    // Optionally carry rho through to the boundary planes as well
    // (UPSTREAM_ISSUES #14: the relaxation constrains only u,v,theta,qv, so
    // the imposed winds carry a mass-flux divergence the local density was
    // never in balance with, and the band absorbs it as spurious w). The
    // frames already hold rho and FillForecastStateMultiFabs has already
    // interpolated it onto the ERF grid -- it was simply dropped here.
    // See HindcastBdyVars.
    bool l_bdy_rho = false, l_mass_consistent = false;
    { ParmParse pp("erf");
      pp.query("hindcast_bdy_rho", l_bdy_rho);
      pp.query("hindcast_mass_consistent_bdy", l_mass_consistent); }
    const int BdyEnd = (l_bdy_rho || l_mass_consistent) ? HindcastBdyVars::NumTypes : BdyEnd0;

    // Same arena convention as init_from_metgrid: CPU+GPU accessible
    Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = The_Pinned_Arena();
#endif

    // Build the plane boxes exactly as init_from_metgrid does
    const auto& lo = geom[lev].Domain().loVect();
    const auto& hi = geom[lev].Domain().hiVect();
    IntVect plo(lo), phi(hi);

    plo[0] = lo[0];              plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = lo[0]+real_width-1; phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xlo(plo, phi);
    Box xlo_plane_no_stag(pbx_xlo);
    Box xlo_plane_x_stag = pbx_xlo; xlo_plane_x_stag.shiftHalf(0,-1);
    Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});

    plo[0] = hi[0]-real_width+1; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0];              phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xhi(plo, phi);
    Box xhi_plane_no_stag(pbx_xhi);
    Box xhi_plane_x_stag = pbx_xhi; xhi_plane_x_stag.shiftHalf(0,1);
    Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});

    plo[1] = lo[1];              plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = lo[1]+real_width-1; phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_ylo(plo, phi);
    Box ylo_plane_no_stag(pbx_ylo);
    Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
    Box ylo_plane_y_stag = pbx_ylo; ylo_plane_y_stag.shiftHalf(1,-1);

    plo[1] = hi[1]-real_width+1; plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = hi[1];              phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_yhi(plo, phi);
    Box yhi_plane_no_stag(pbx_yhi);
    Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
    Box yhi_plane_y_stag = pbx_yhi; yhi_plane_y_stag.shiftHalf(1,1);

    bdy_data_xlo.resize(ntimes);
    bdy_data_xhi.resize(ntimes);
    bdy_data_ylo.resize(ntimes);
    bdy_data_yhi.resize(ntimes);

    for (int itime(0); itime < ntimes; itime++) {
        bdy_data_xlo[itime].resize(BdyEnd);
        bdy_data_xhi[itime].resize(BdyEnd);
        bdy_data_ylo[itime].resize(BdyEnd);
        bdy_data_yhi[itime].resize(BdyEnd);
        for (int nvar(0); nvar<BdyEnd; ++nvar) {
            if (nvar==MetGridBdyVars::U) {
                bdy_data_xlo[itime][nvar].resize(xlo_plane_x_stag, 1, Arena_Used);
                bdy_data_xhi[itime][nvar].resize(xhi_plane_x_stag, 1, Arena_Used);
                bdy_data_ylo[itime][nvar].resize(ylo_plane_x_stag, 1, Arena_Used);
                bdy_data_yhi[itime][nvar].resize(yhi_plane_x_stag, 1, Arena_Used);
            } else if (nvar==MetGridBdyVars::V) {
                bdy_data_xlo[itime][nvar].resize(xlo_plane_y_stag, 1, Arena_Used);
                bdy_data_xhi[itime][nvar].resize(xhi_plane_y_stag, 1, Arena_Used);
                bdy_data_ylo[itime][nvar].resize(ylo_plane_y_stag, 1, Arena_Used);
                bdy_data_yhi[itime][nvar].resize(yhi_plane_y_stag, 1, Arena_Used);
            } else {
                bdy_data_xlo[itime][nvar].resize(xlo_plane_no_stag, 1, Arena_Used);
                bdy_data_xhi[itime][nvar].resize(xhi_plane_no_stag, 1, Arena_Used);
                bdy_data_ylo[itime][nvar].resize(ylo_plane_no_stag, 1, Arena_Used);
                bdy_data_yhi[itime][nvar].resize(yhi_plane_no_stag, 1, Arena_Used);
            }
        }

        // Interpolate this frame onto the ERF grid (forecast_state_1 is
        // scratch here; the run-time machinery re-reads its frames at the
        // first Evolve step).
        FillForecastStateMultiFabs(lev, bin_files[itime], z_phys_nd[lev], forecast_state_1);

        const MultiFab& fcons = forecast_state_1[lev][Vars::cons];
        const MultiFab& fxvel = forecast_state_1[lev][Vars::xvel];
        const MultiFab& fyvel = forecast_state_1[lev][Vars::yvel];

        for (int nvar(0); nvar<BdyEnd; ++nvar) {
            const MultiFab* src = nullptr;
            int scomp = 0;
            if      (nvar==MetGridBdyVars::U)  { src = &fxvel; scomp = 0; }
            else if (nvar==MetGridBdyVars::V)  { src = &fyvel; scomp = 0; }
            else if (nvar==MetGridBdyVars::T)  { src = &fcons; scomp = RhoTheta_comp; } // plain theta
            else if (nvar==MetGridBdyVars::QV) { src = &fcons; scomp = RhoQ1_comp;    } // plain qv
            else if (nvar==HindcastBdyVars::RHO) { src = &fcons; scomp = Rho_comp;    } // density
            strip_to_global_fab(*src, scomp, bdy_data_xlo[itime][nvar].box(), bdy_data_xlo[itime][nvar]);
            strip_to_global_fab(*src, scomp, bdy_data_xhi[itime][nvar].box(), bdy_data_xhi[itime][nvar]);
            strip_to_global_fab(*src, scomp, bdy_data_ylo[itime][nvar].box(), bdy_data_ylo[itime][nvar]);
            strip_to_global_fab(*src, scomp, bdy_data_yhi[itime][nvar].box(), bdy_data_yhi[itime][nvar]);
        }
    } // itime

    Print() << "HindCast real BCs: filled " << ntimes << " boundary-plane times "
            << "(width " << real_width << " cells, interval "
            << bdy_time_interval << " s) from " << folder << std::endl;
}
#endif

/**
 * Rebuild the HSE base state and the thermodynamic state from the
 * time-interpolated ERA5 frame at initialization time.
 *
 * Without this, the HindCast pathway starts from the problem default -- a dry
 * 300 K isentrope -- and the theta/qv lateral sponge then relaxes the bands
 * toward ERA5 theta (330-430 K aloft) against a 300 K interior/base state:
 * measured result is O(30-100 K) standing lateral theta contrasts aloft,
 * w of +-50 m/s within ~15 model minutes, and NaN. The interior and base
 * state must START on the observed stratification.
 *
 * Method: the interpolated ERA5 density becomes the base-state density;
 * erf_enforce_hse integrates the hydrostatic pressure from it and derives
 * the consistent theta -- so base and state agree exactly (zero buoyancy at
 * init) and carry the real stratification. qv comes from the frame. Momenta
 * remain zero (spin-up from rest, as before).
 */
void
ERF::init_thermo_from_hindcast (const int lev)
{
    // The forecast state was first filled during init_stuff, BEFORE the
    // terrain arrays were built (z_phys_nd was still zero). Re-run the
    // interpolation with the final grid; regrid_forces_file_read=true forces
    // the re-read without advancing the frame clock.
    WeatherDataInterpolation(lev, t_new[lev], z_phys_nd, true);

    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp , 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp , 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    MultiFab& fcons = forecast_state_interp[lev][Vars::cons];
    MultiFab& cons  = vars_new[lev][Vars::cons];

    const bool l_has_moist = (solverChoice.moisture_type != MoistureType::None) &&
                             (cons.nComp() > RhoQ1_comp) && (fcons.nComp() > RhoQ1_comp);

    IntVect ngv = r_hse.nGrowVect();
    ngv.min(fcons.nGrowVect());
    MultiFab::Copy(r_hse, fcons, Rho_comp, 0, 1, ngv);
    if (l_has_moist) {
        // Forecast state stores PLAIN qv (velocity convention)
        MultiFab::Copy(qv_hse, fcons, RhoQ1_comp, 0, 1, ngv);
    }

    erf_enforce_hse(lev, r_hse, p_hse, pi_hse, th_hse, qv_hse, z_phys_cc[lev]);
    (*physbcs_base[lev])(base_state[lev], 0, base_state[lev].nComp(), base_state[lev].nGrowVect());

    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const Box& gbx = mfi.growntilebox(1);
        const Array4<Real      >& cons_arr = cons.array(mfi);
        const Array4<Real const>& r_arr    = r_hse.const_array(mfi);
        const Array4<Real const>& th_arr   = th_hse.const_array(mfi);
        const Array4<Real const>& f_arr    = fcons.const_array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cons_arr(i,j,k,Rho_comp)      = r_arr(i,j,k);
            cons_arr(i,j,k,RhoTheta_comp) = r_arr(i,j,k) * th_arr(i,j,k);
            if (l_has_moist) {
                cons_arr(i,j,k,RhoQ1_comp) = r_arr(i,j,k) * f_arr(i,j,k,RhoQ1_comp);
            }
        });
    }

    Print() << "HindCast init: base state and thermodynamic state rebuilt from "
            << "the interpolated ERA5 frame (theta/qv coupling); lev " << lev
            << " rho min/max " << cons.min(Rho_comp) << " " << cons.max(Rho_comp) << std::endl;

    hindcast_check_mass_consistency(lev);
}

/**
 * Measure the discrete continuity residual of the ERA5 target field under ERF's
 * OWN operator, i.e. the one in AdvectionSrcForRho:
 *
 *     (detJ/m^2) drho/dt + dFx/dx + dFy/dy + dFz/dzeta = 0,
 *     Fx = ax*(rho u)/mf_uy,  Fy = ay*(rho v)/mf_vx,  Fz = az*Omega/m^2.
 *
 * Omega is hard-zeroed at k = klo and the lid is a SlipWall, so Fz telescopes
 * out of a column sum and
 *
 *     R(i,j) = sum_k [ (detJ/m^2) d(rho_t)/dt + d(Fx_t)/dx + d(Fy_t)/dy ] dzeta
 *
 * (subscript t = the ERA5 target field) must vanish for ANY field the band can
 * be relaxed toward -- it is the compatibility condition of the two-point
 * problem for Omega_t. Integrating upward from Omega_t = 0 gives the vertical
 * mass flux the imposed horizontal fluxes demand; |Omega_t|/rho_t is its
 * equivalent vertical velocity, directly
 * comparable to the w > 1 m/s band criterion. R != 0 means no relaxation
 * target can be mass-consistent and the target divergence itself needs the
 * column-mass correction.
 */
void
ERF::hindcast_check_mass_consistency (const int lev)
{
    const Real dT = solverChoice.hindcast_data_interval_in_hrs * Real(3600.0);
    if (dT <= zero) { return; }

    const MultiFab& fc1  = forecast_state_1[lev][Vars::cons];
    const MultiFab& fc2  = forecast_state_2[lev][Vars::cons];
    const MultiFab& ftgt = forecast_state_interp[lev][Vars::cons];
    const MultiFab& xvel = vars_new[lev][Vars::xvel];
    const MultiFab& yvel = vars_new[lev][Vars::yvel];

    const auto dxInv  = geom[lev].InvCellSizeArray();
    const Real dzeta  = geom[lev].CellSize(2);
    const Box& domain = geom[lev].Domain();
    const auto dom_lo = lbound(domain);
    const auto dom_hi = ubound(domain);
    const int  bw     = (real_width > 0) ? real_width : 1;

    // (band |Omega*|/rho* max, band top-residual max, interior same two)
    ReduceOps<ReduceOpMax,ReduceOpMax,ReduceOpMax,ReduceOpMax> reduce_op;
    ReduceData<Real,Real,Real,Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    bool z_split = false;
    for (MFIter mfi(ftgt); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        // The upward recursion needs the whole column on one box.
        if (bx.smallEnd(2) != dom_lo.z || bx.bigEnd(2) != dom_hi.z) { z_split = true; continue; }
        const Box bx2 = makeSlab(bx,2,bx.smallEnd(2));

        const Array4<const Real>& r1  = fc1.const_array(mfi);
        const Array4<const Real>& r2  = fc2.const_array(mfi);
        const Array4<const Real>& rt  = ftgt.const_array(mfi);
        const Array4<const Real>& u   = xvel.const_array(mfi);
        const Array4<const Real>& v   = yvel.const_array(mfi);
        const Array4<const Real>& axa = ax[lev]->const_array(mfi);
        const Array4<const Real>& aya = ay[lev]->const_array(mfi);
        const Array4<const Real>& dJ  = detJ_cc[lev]->const_array(mfi);
        const Array4<const Real>& mfx = mapfac[lev][MapFacType::m_x]->const_array(mfi);
        const Array4<const Real>& mfy = mapfac[lev][MapFacType::m_y]->const_array(mfi);
        const Array4<const Real>& mfu = mapfac[lev][MapFacType::u_y]->const_array(mfi);
        const Array4<const Real>& mfv = mapfac[lev][MapFacType::v_x]->const_array(mfi);

        const int klo = dom_lo.z, khi = dom_hi.z;

        reduce_op.eval(bx2, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> ReduceTuple
        {
            // Clamp into the valid region: the target's physical-boundary
            // ghosts are zero-initialized, and the wall face uses the
            // zero-gradient density that fill_from_realbdy imposes.
            auto rho_at = [=] (int ii, int jj, int kk) {
                ii = amrex::min(amrex::max(ii,dom_lo.x),dom_hi.x);
                jj = amrex::min(amrex::max(jj,dom_lo.y),dom_hi.y);
                return rt(ii,jj,kk,Rho_comp);
            };

            const Real msq = mfx(i,j,0) * mfy(i,j,0);
            Real Fz   = Real(0.0);   // az*Omega_t/m^2, vanishes at the terrain
            Real wmax = Real(0.0);

            for (int k = klo; k <= khi; ++k) {
                const Real rc = rho_at(i,j,k);

                const Real rux_lo = u(i  ,j,k) * myhalf*(rho_at(i-1,j,k) + rc);
                const Real rux_hi = u(i+1,j,k) * myhalf*(rc + rho_at(i+1,j,k));
                const Real rvy_lo = v(i,j  ,k) * myhalf*(rho_at(i,j-1,k) + rc);
                const Real rvy_hi = v(i,j+1,k) * myhalf*(rc + rho_at(i,j+1,k));

                const Real Fx_lo = axa(i  ,j,k) * rux_lo / mfu(i  ,j,0);
                const Real Fx_hi = axa(i+1,j,k) * rux_hi / mfu(i+1,j,0);
                const Real Fy_lo = aya(i,j  ,k) * rvy_lo / mfv(i,j  ,0);
                const Real Fy_hi = aya(i,j+1,k) * rvy_hi / mfv(i,j+1,0);

                const Real drdt = (r2(i,j,k,Rho_comp) - r1(i,j,k,Rho_comp)) / dT;

                const Real res = dJ(i,j,k)/msq * drdt
                               + (Fx_hi - Fx_lo) * dxInv[0]
                               + (Fy_hi - Fy_lo) * dxInv[1];

                Fz -= res * dzeta;               // Fz_{k+1} = Fz_k - dzeta*res
                const Real w_eq = (Fz * msq) / amrex::max(rc, Real(1.e-6));
                wmax = amrex::max(wmax, std::abs(w_eq));
            }

            // Fz after the loop is the flux the lid would have to pass; the
            // SlipWall forbids it, so this is the incompatibility.
            const Real w_top = std::abs((Fz * msq) / amrex::max(rho_at(i,j,khi), Real(1.e-6)));

            const bool in_band = (i < dom_lo.x + bw) || (i > dom_hi.x - bw) ||
                                 (j < dom_lo.y + bw) || (j > dom_hi.y - bw);

            return { in_band ? wmax : Real(0.0), in_band ? w_top : Real(0.0),
                     in_band ? Real(0.0) : wmax, in_band ? Real(0.0) : w_top };
        });
    }

    ReduceTuple hv = reduce_data.value(reduce_op);
    Real band_w = amrex::get<0>(hv), band_top = amrex::get<1>(hv);
    Real int_w  = amrex::get<2>(hv), int_top  = amrex::get<3>(hv);
    ParallelDescriptor::ReduceRealMax(band_w);
    ParallelDescriptor::ReduceRealMax(band_top);
    ParallelDescriptor::ReduceRealMax(int_w);
    ParallelDescriptor::ReduceRealMax(int_top);

    if (z_split) {
        Warning("hindcast mass-consistency check skipped boxes split in z");
    }

    // The model's rho is the HSE-rebalanced density (erf_enforce_hse above),
    // not the raw interpolated ERA5 density. Relaxing rho toward the raw field
    // would drive the state off its own discrete hydrostatic balance, so
    // report the discrepancy: it is the difference between "mass-consistent"
    // and "hydrostatically consistent" targets.
    Real rdiff_max = zero, rdiff_l1 = zero, rsum = zero;
    {
        MultiFab diff(ftgt.boxArray(), ftgt.DistributionMap(), 1, 0);
        MultiFab::Copy   (diff, ftgt, Rho_comp, 0, 1, 0);
        MultiFab::Subtract(diff, vars_new[lev][Vars::cons], Rho_comp, 0, 1, 0);
        rdiff_max = diff.norm0();
        rdiff_l1  = diff.norm1();
        rsum      = vars_new[lev][Vars::cons].norm1(Rho_comp);
    }
    Print() << "[mass-consistency] target rho vs model (HSE-rebalanced) rho: max|drho| = "
            << rdiff_max << " kg/m^3, mean|drho|/mean(rho) = "
            << ((rsum > zero) ? rdiff_l1/rsum : zero) << std::endl;
    Print() << "[mass-consistency] ERA5 target under ERF's discrete continuity operator, lev "
            << lev << ":\n"
            << "    band (" << bw << " cells): max implied |w| = " << band_w
            << " m/s, max lid incompatibility = " << band_top << " m/s\n"
            << "    interior          : max implied |w| = " << int_w
            << " m/s, max lid incompatibility = " << int_top << " m/s" << std::endl;
}

void
ERF::WeatherDataInterpolation(const int lev,
                              const Real time,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd,
                              bool regrid_forces_file_read)
{

    static amrex::Vector<Real> next_read_forecast_time;
    static amrex::Vector<Real> last_read_forecast_time;

    const int nlevs = static_cast<int>(a_z_phys_nd.size());

    Real hindcast_data_interval = solverChoice.hindcast_data_interval_in_hrs*Real(3600.0);

    // Initialize static vectors once
    if (next_read_forecast_time.empty()) {
        next_read_forecast_time.resize(nlevs, -one);
        last_read_forecast_time.resize(nlevs, -one);
        Print() << "Initializing the time vector values here by " << lev << std::endl;
    }

    if (next_read_forecast_time[lev] < zero) {
        int next_multiple = static_cast<int>(time / hindcast_data_interval);
        next_read_forecast_time[lev] = next_multiple * hindcast_data_interval;
        last_read_forecast_time[lev] = next_read_forecast_time[lev];
    }

    if (time >= next_read_forecast_time[lev] or regrid_forces_file_read) {

        Print() << "Data reading happening at level " << lev << std::endl;

        std::string folder = solverChoice.hindcast_boundary_data_dir;

        // Check if folder exists and is a directory
        if (!fs::exists(folder) || !fs::is_directory(folder)) {
            throw std::runtime_error("Error: Folder '" + folder + "' does not exist or is not a directory.");
        }

        std::vector<std::string> bin_files;

        for (const auto& entry : fs::directory_iterator(folder)) {
            if (!entry.is_regular_file()) continue;

            std::string fname = entry.path().filename().string();
            if (fname.size() >= 4 && fname.substr(fname.size() - 4) == ".bin") {
                bin_files.push_back(entry.path().string());
            }
        }
        std::sort(bin_files.begin(), bin_files.end());

    // Check if no .bin files were found
        if (bin_files.empty()) {
            throw std::runtime_error("Error: No .bin files found in folder '" + folder + "'.");
        }

        std::string filename1, filename2;

        int idx1 = static_cast<int>(time / hindcast_data_interval);
        int idx2 = static_cast<int>(time / hindcast_data_interval)+1;
        Print() << "Reading weather data " << time << " " << idx1 << " " << idx2 <<" " << bin_files.size() << std::endl;

        if (idx2 >= static_cast<int>(bin_files.size())) {
            throw std::runtime_error("Error: Not enough .bin files to cover time " + std::to_string(time));
        }

        filename1 = bin_files[idx1];
        filename2 = bin_files[idx2];

        FillForecastStateMultiFabs(lev, filename1, a_z_phys_nd[lev], forecast_state_1);
        FillForecastStateMultiFabs(lev, filename2, a_z_phys_nd[lev], forecast_state_2);

         // Create the time-interpolated forecast state
        //CreateForecastStateMultiFabs(forecast_state_interp);
        if(!regrid_forces_file_read){
            last_read_forecast_time[lev] = next_read_forecast_time[lev];
            next_read_forecast_time[lev] += hindcast_data_interval;
            Print() << "Next forecast time getting updated here " << std::endl;
        }
    }

    Real prev_read_time = last_read_forecast_time[lev];
    Real alpha1 = one - (time - prev_read_time)/hindcast_data_interval;
    Real alpha2 = one - alpha1;

    amrex::Print()<< "The values of alpha1 and alpha2 are " << alpha1 << " "<< alpha2 <<std::endl;

    if (alpha1 < zero || alpha1 > one ||
    alpha2 < zero || alpha2 > one)
    {
        std::stringstream ss;
        ss << "Interpolation weights for hindcast files are incorrect: "
        << "alpha1 = " << alpha1 << ", alpha2 = " << alpha2;
       Abort(ss.str());
    }

    MultiFab& erf_mf_cons   = forecast_state_interp[lev][Vars::cons];
    MultiFab& erf_mf_xvel   = forecast_state_interp[lev][Vars::xvel];
    MultiFab& erf_mf_yvel   = forecast_state_interp[lev][Vars::yvel];
    //MultiFab& erf_mf_zvel   = forecast_state_interp[0][Vars::zvel];
    MultiFab& erf_mf_latlon = forecast_state_interp[lev][4];

    // Fill the time-interpolated forecast states
    MultiFab::LinComb(forecast_state_interp[lev][Vars::cons],
                      alpha1, forecast_state_1[lev][Vars::cons], 0,
                      alpha2, forecast_state_2[lev][Vars::cons], 0,
                      0, erf_mf_cons.nComp(), forecast_state_interp[lev][Vars::cons].nGrow());
    // The velocity fill kernels above cover VALID face boxes only (no ghosts),
    // so forecast_state_1/2 velocity ghost cells are uninitialized. Blend over
    // 0 ghosts and fill the interp ghosts by exchange instead: blending over
    // nGrow copies uninitialized memory into the interp state, which the
    // hindcast sponge kernels then read at box perimeters (NaN in
    // memory-dependent regions).
    MultiFab::LinComb(forecast_state_interp[lev][Vars::xvel],
                      alpha1, forecast_state_1[lev][Vars::xvel], 0,
                      alpha2, forecast_state_2[lev][Vars::xvel], 0,
                      0, erf_mf_xvel.nComp(), 0);
    MultiFab::LinComb(forecast_state_interp[lev][Vars::yvel],
                      alpha1, forecast_state_1[lev][Vars::yvel], 0,
                      alpha2, forecast_state_2[lev][Vars::yvel], 0,
                      0, erf_mf_yvel.nComp(), 0);
    forecast_state_interp[lev][Vars::xvel].FillBoundary(geom[lev].periodicity());
    forecast_state_interp[lev][Vars::yvel].FillBoundary(geom[lev].periodicity());
    MultiFab::LinComb(forecast_state_interp[lev][4],
                      alpha1, forecast_state_1[lev][4], 0,
                      alpha2, forecast_state_2[lev][4], 0,
                      0, erf_mf_latlon.nComp(), forecast_state_interp[lev][4].nGrow());

    /*Vector<std::string> varnames_plot_mf = {
    "rho", "rhotheta", "rhoqv", "rhoqc", "rhoqr", "xvel", "yvel", "zvel", "latitude", "longitude"
    }; // Customize variable names

    std::string pltname = "plt_interp";

    MultiFab plot_mf(erf_mf_cons.boxArray(), erf_mf_cons.DistributionMap(),
                     10, 0);

    plot_mf.setVal(0.0);

    for (MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
        const Array4<Real> &plot_mf_arr = plot_mf.array(mfi);
        const Array4<Real> &erf_mf_cons_arr = erf_mf_cons.array(mfi);
        const Array4<Real> &erf_mf_xvel_arr = erf_mf_xvel.array(mfi);
        const Array4<Real> &erf_mf_yvel_arr = erf_mf_yvel.array(mfi);
        const Array4<Real> &erf_mf_zvel_arr = erf_mf_zvel.array(mfi);
        const Array4<Real> &erf_mf_latlon_arr = erf_mf_latlon.array(mfi);

        const Box& bx = mfi.validbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plot_mf_arr(i,j,k,0) = erf_mf_cons_arr(i,j,k,Rho_comp);
            plot_mf_arr(i,j,k,1) = erf_mf_cons_arr(i,j,k,RhoTheta_comp);
            plot_mf_arr(i,j,k,2) = erf_mf_cons_arr(i,j,k,RhoQ1_comp);
            plot_mf_arr(i,j,k,3) = erf_mf_cons_arr(i,j,k,RhoQ2_comp);
            plot_mf_arr(i,j,k,4) = erf_mf_cons_arr(i,j,k,RhoQ3_comp);

            plot_mf_arr(i,j,k,5) = (erf_mf_xvel_arr(i,j,k,0) + erf_mf_xvel_arr(i+1,j,k,0))/two;
            plot_mf_arr(i,j,k,6) = (erf_mf_yvel_arr(i,j,k,0) + erf_mf_yvel_arr(i,j+1,k,0))/two;
            plot_mf_arr(i,j,k,7) = (erf_mf_zvel_arr(i,j,k,0) + erf_mf_zvel_arr(i,j,k+1,0))/two;

            plot_mf_arr(i,j,k,8) = erf_mf_latlon_arr(i,j,k,0);
            plot_mf_arr(i,j,k,9) = erf_mf_latlon_arr(i,j,k,1);
        });
    }


    WriteSingleLevelPlotfile(
            pltname,
            plot_mf,
            varnames_plot_mf,
            geom[0],
            time,
            0 // level
        );*/
}

#endif
