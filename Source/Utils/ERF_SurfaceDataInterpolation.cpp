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

void
ERF::FillSurfaceStateMultiFabs(const int lev,
                               const std::string& filename,
                               Vector<MultiFab>& surface_state)
{
       // Open the binary file in input mode
    std::ifstream infile(filename, std::ios::binary);
    if (!infile) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
    }
    Vector<Real> xvec_h, yvec_h, zvec_h;
    Vector<Real> sst_h, q_star_h, t_star_h, u_star_h, ls_mask_h;
    Vector<Real> alb_h, junk_h;   // field 5 = forecast albedo (fal patch); junk guards unknown extras

    int nx, ny, nz, ndata;
    float value;

    // Read the four integers
    infile.read(reinterpret_cast<char*>(&nx), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ny), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nz), sizeof(int));
    infile.read(reinterpret_cast<char*>(&ndata), sizeof(int));

    amrex::Gpu::DeviceVector<Real> xvec_d(nx*ny*nz), yvec_d(nx*ny*nz), zvec_d(nz);

    // The frame format carries two nx*ny lat/lon planes between the header and
    // the x/y/z coordinate vectors. The 3D path reads them (ReadCustomBinaryIC
    // takes latvec_h/lonvec_h); this parser did not skip them, so every read
    // below was offset by 2*nx*ny - (nx+ny+nz) floats: xvec_h held latitudes,
    // sst_h held longitudes, and alb_h ran past EOF. Targets then clamped to a
    // single corner in bilinear_interpolation_2d and no SST stencil point ever
    // passed the 271-305 K test, so the whole surface state came back as
    // sentinels (mask 0, SST -1, albedo at its 0.03 floor).
    infile.seekg(std::streamoff(2)*nx*ny*std::streamoff(sizeof(float)), std::ios::cur);

    for(int i=0; i<nx; i++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        xvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, xvec_h.begin(), xvec_h.end(), xvec_d.begin());

    for(int j=0; j<ny; j++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        yvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, yvec_h.begin(), yvec_h.end(), yvec_d.begin());

    for(int k=0; k<nz; k++) {
        infile.read(reinterpret_cast<char*>(&value), sizeof(float));
        zvec_h.emplace_back(value);
    }
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, zvec_h.begin(), zvec_h.end(), zvec_d.begin());

         // Vector to store the data

    Vector<Real>* data_h = nullptr; // Declare pointer outside the loop

    Real* xvec_d_ptr = xvec_d.data();
    Real* yvec_d_ptr = yvec_d.data();

    Real dxvec = (xvec_h[nx-1]-xvec_h[0])/(nx-1);
    Real dyvec = (yvec_h[ny-1]-yvec_h[0])/(ny-1);

    // Read the file
    for(int idx=0; idx<ndata; idx++){
        if(idx == 0){
            data_h = &sst_h;
        } else if (idx==1) {
            data_h = &q_star_h;
        } else if (idx==2) {
            data_h = &t_star_h;
        } else if (idx==3) {
            data_h = &u_star_h;
        } else if(idx==4) {
            data_h = &ls_mask_h;
        } else if(idx==5) {
            data_h = &alb_h;
        } else {
            // Unknown extra fields must not fall through into the last
            // named vector (that silently corrupts it) -- park them.
            junk_h.clear();
            data_h = &junk_h;
        }
        for(int k=0; k<nz; k++) {
            for(int j=0; j<ny; j++) {
                for(int i=0; i<nx; i++) {
                    infile.read(reinterpret_cast<char*>(&value), sizeof(float));
                    //if(idx == 3) {
                        //printf("theta is %0.15g, %0.15g, %0.15g %0.15g\n", xvec_h[i], yvec_h[j], zvec_h[k], value);
                    //}
                    data_h->emplace_back(value);
                }
            }
        }
    }

    infile.close();

    // Audit the SOURCE array before it is interpolated, so we can separate what the
    // generator wrote from what our own interpolation creates.
    {
        int n_water = 0, n_bad_water = 0, n_bad_land = 0;
        for (size_t m = 0; m < sst_h.size(); ++m) {
            const bool water = (m < ls_mask_h.size()) && (ls_mask_h[m] < Real(0.5));
            const bool bad   = !(sst_h[m] > Real(271.0) && sst_h[m] < Real(305.0));
            if (water) { ++n_water; if (bad) ++n_bad_water; }
            else if (bad) { ++n_bad_land; }
        }
        static bool src_reported = false;
        if (!src_reported) {
            src_reported = true;
            Print() << "HindCast SST source frame: " << n_water << " water points, "
                    << n_bad_water << " of them outside 271-305 K (generator smear); "
                    << n_bad_land << " land points carrying fills" << std::endl;
        }
    }

    amrex::Gpu::DeviceVector<Real> ls_mask_d(nx*ny*nz), sst_d(nx*ny*nz);

    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, ls_mask_h.begin(), ls_mask_h.end(), ls_mask_d.begin());
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, sst_h.begin(), sst_h.end(), sst_d.begin());

    // Forecast albedo (field 5, fal patch): absent in pre-patch 5-field
    // frames -- signal that with an empty device vector.
    const bool have_alb = (static_cast<int>(alb_h.size()) == nx*ny*nz);
    amrex::Gpu::DeviceVector<Real> alb_d(have_alb ? nx*ny*nz : 0);
    if (have_alb) {
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, alb_h.begin(), alb_h.end(), alb_d.begin());
    }

    Real* ls_mask_d_ptr = ls_mask_d.data();
    Real* sst_d_ptr   = sst_d.data();
    Real* alb_d_ptr   = have_alb ? alb_d.data() : nullptr;

    // Use the GpuArray accessors: ProbLo()/CellSize() return host pointers,
    // which are invalid when captured into a device lambda (CUDA error 700).
    // The 3D path (ERF_WeatherDataInterpolation.cpp) already does this.
    const auto prob_lo  = geom[lev].ProbLoArray();
    const auto dx       = geom[lev].CellSizeArray();

    for (amrex::MFIter mfi(surface_state[lev]); mfi.isValid(); ++mfi) {
        const Box gbx = mfi.growntilebox();
        const Array4<Real>& surf_arr = surface_state[lev].array(mfi);

        ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

            // NOTE: write EVERY k-plane of the slab, not just k==0. The 2D
            // surface MultiFab inherits the 3D state's ghost vector, so it
            // owns z-ghost planes (k != 0); leaving them unwritten fed
            // uninitialized memory into every downstream copy/blend of this
            // field (flagged by compute-sanitizer initcheck). The field is
            // k-independent, so all planes take the same value.
            const Real x        = prob_lo[0] + (i + myhalf) * dx[0];
            const Real y        = prob_lo[1] + (j + myhalf) * dx[1];

            // First interpolate where the weather data is available from
            Real tmp_ls_mask, tmp_sst;

            bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                      dxvec, dyvec,
                                      nx, ny,
                                      x, y,
                                      ls_mask_d_ptr, tmp_ls_mask);

            // MASKED bilinear for SST. The raw call blends across the coastline into
            // 9999 land fills, and on the finer ERF mesh that produces a continuum of
            // intermediate values -- including ones that look physical (measured:
            // 62.4 C passing a 200-340 K guard). Accumulate only source points that
            // are BOTH water AND in physical range, renormalising the weights, so no
            // contaminated value is ever constructed. A land-masked point can still
            // carry a fill, and a water point can still be smeared in the source, so
            // both criteria are required.
            //
            // Cells whose stencil contains no valid source are flagged (-1) rather
            // than given a fabricated value; the nearest-valid fill downstream is what
            // handles those, which is what it is for.
            {
                const Real fi = (x - xvec_d_ptr[0]) / dxvec;
                const Real fj = (y - yvec_d_ptr[0]) / dyvec;
                int i0 = static_cast<int>(std::floor(fi));
                int j0 = static_cast<int>(std::floor(fj));
                const Real wx = fi - static_cast<Real>(i0);
                const Real wy = fj - static_cast<Real>(j0);
                Real num = Real(0.0), den = Real(0.0);
                for (int dj = 0; dj < 2; ++dj) {
                for (int di = 0; di < 2; ++di) {
                    const int ii = amrex::min(amrex::max(i0+di, 0), nx-1);
                    const int jj = amrex::min(amrex::max(j0+dj, 0), ny-1);
                    const Real sv = sst_d_ptr    [jj*nx + ii];
                    const Real lv = ls_mask_d_ptr[jj*nx + ii];
                    if (lv < Real(0.5) && sv > Real(271.0) && sv < Real(305.0)) {
                        const Real w = (di ? wx : Real(1.0)-wx) * (dj ? wy : Real(1.0)-wy);
                        num += w*sv; den += w;
                    }
                }}
                tmp_sst = (den > Real(1.0e-8)) ? num/den : Real(-1.0);
            }

            surf_arr(i, j, k, 0) = std::min(tmp_ls_mask, amrex::Real(1.0));
            surf_arr(i, j, k, 1) = tmp_sst;

            // comp 3 = LAND skin temperature: the mirror of the stencil above.
            // The frames carry TSK over land in this same field (write_frame
            // stores sst = where(water & valid, SST, TSK)), but the water-only
            // mask above throws it away, so it has to be interpolated a second
            // time with the opposite mask. Without this, land cells of comp 1
            // are the -1 no-valid-source flag, not a temperature.
            //
            // The physical band is wide (200-340 K), not the ocean's 271-305:
            // real land skin temperature on this domain reaches 261 K on the
            // high terrain, and an ocean band would reject exactly the cells
            // that carry the signal. -1 means "no land source point in the
            // stencil", handled downstream, never a fabricated value.
            {
                const Real fi = (x - xvec_d_ptr[0]) / dxvec;
                const Real fj = (y - yvec_d_ptr[0]) / dyvec;
                int i0 = static_cast<int>(std::floor(fi));
                int j0 = static_cast<int>(std::floor(fj));
                const Real wx = fi - static_cast<Real>(i0);
                const Real wy = fj - static_cast<Real>(j0);
                Real num = Real(0.0), den = Real(0.0);
                for (int dj = 0; dj < 2; ++dj) {
                for (int di = 0; di < 2; ++di) {
                    const int ii = amrex::min(amrex::max(i0+di, 0), nx-1);
                    const int jj = amrex::min(amrex::max(j0+dj, 0), ny-1);
                    const Real sv = sst_d_ptr    [jj*nx + ii];
                    const Real lv = ls_mask_d_ptr[jj*nx + ii];
                    if (lv >= Real(0.5) && sv > Real(200.0) && sv < Real(340.0)) {
                        const Real w = (di ? wx : Real(1.0)-wx) * (dj ? wy : Real(1.0)-wy);
                        num += w*sv; den += w;
                    }
                }}
                surf_arr(i, j, k, 3) = (den > Real(1.0e-8)) ? num/den : Real(-1.0);
            }

            // comp 2 = surface albedo (fal); -1 marks "not in this file"
            // so downstream wiring can fall back to constants.
            if (alb_d_ptr) {
                Real tmp_alb;
                bilinear_interpolation_2d(xvec_d_ptr, yvec_d_ptr,
                                          dxvec, dyvec,
                                          nx, ny,
                                          x, y,
                                          alb_d_ptr, tmp_alb);
                // Clamp to a physical range (bilinear blending across the
                // coastline and GRIB packing can nick the edges).
                surf_arr(i, j, k, 2) = amrex::min(amrex::max(tmp_alb, Real(0.03)), Real(0.95));
            } else {
                surf_arr(i, j, k, 2) = Real(-1.0);
            }
        });
    }

    // LIFETIME BARRIER -- same defect as FillForecastStateMultiFabs, see
    // UPSTREAM_ISSUES 29. xvec_d/yvec_d/zvec_d/ls_mask_d/sst_d/alb_d are
    // function-local DeviceVectors whose raw pointers are captured by the
    // ParallelFors above; returning without a barrier frees them into the arena
    // while those kernels may still be reading.
    amrex::Gpu::streamSynchronize();
}

void
ERF::SurfaceDataInterpolation(const int lev,
                              const Real time,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd,
                              bool regrid_forces_file_read)
{

    static amrex::Vector<Real> next_read_forecast_time;
    static amrex::Vector<Real> last_read_forecast_time;

    const int nlevs = a_z_phys_nd.size();

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

        std::string folder = solverChoice.hindcast_surface_data_dir;

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
        Print() << "Reading surface data " << time << " " << idx1 << " " << idx2 <<" " << bin_files.size() << std::endl;

        if (idx2 >= static_cast<int>(bin_files.size())) {
            throw std::runtime_error("Error: Not enough .bin files to cover time " + std::to_string(time));
        }

        filename1 = bin_files[idx1];
        filename2 = bin_files[idx2];

        FillSurfaceStateMultiFabs(lev, filename1, surface_state_1);
        FillSurfaceStateMultiFabs(lev, filename2, surface_state_2);

        // Wire the hindcast surface data into the MOST surface layer:
        // an LSM (or MYNN-EDMF's surface layer) over a domain with water
        // requires SST (m_sst_lev) and a land mask; previously only the
        // WRF lower-boundary path populated these. Register a single
        // time slot (the surface layer's nt=1 branch uses it directly)
        // and refresh its data at every 3-hourly read. The surface layer
        // stores raw pointers at construction, so the allocation must
        // persist and only its contents change.
        {
            MultiFab& surf_mf = surface_state_1[lev];
            if (sst_lev[lev].empty()) {
                sst_lev[lev].resize(1);
                sst_lev[lev][0] = std::make_unique<MultiFab>(
                    surf_mf.boxArray(), surf_mf.DistributionMap(), 1,
                    surf_mf.nGrowVect());
            }
            // TSK carries the COMPLETE surface temperature field: the land skin
            // temperature AND, over water, the SST. That is upstream's contract,
            // not a local convention -- ERF_ReadFromWRFLow.cpp:238 copies SST
            // into TSK over water for the same reason, and SurfaceLayer sets
            // m_ignore_sst as soon as TSK exists (ERF_SurfaceLayer.H:465) so
            // EVERY cell reads it. conus404_to_bin.py's write_frame already
            // performs that merge offline
            //     sst = np.where(water & ~bad, SST, TSK)
            // so frame field 1 is exactly the field TSK wants.
            //
            // Without this, land t_surf fell through to default_land_surf_temp
            // (or, with an LSM registered, to get_lsm_tsurf's soil constant):
            // one number for the whole run. Measured on the 2020-12-28 frames,
            // the real land skin temperature swings 276.8-285.9 K in the domain
            // mean over a day -- 261-296 K across individual cells -- minimum at
            // dawn, maximum in early afternoon. RRTMGP takes its radiating skin
            // temperature from the same t_surf, so the constant removed the land
            // diurnal cycle from the longwave as well.
            if (tsk_lev[lev].empty()) { tsk_lev[lev].resize(1); }
            if (!tsk_lev[lev][0]) {
                tsk_lev[lev][0] = std::make_unique<MultiFab>(
                    surf_mf.boxArray(), surf_mf.DistributionMap(), 1,
                    surf_mf.nGrowVect());
            }
            // The default all-land lmask is allocated later in init_stuff;
            // at the first (init-time) call here it does not exist yet, so
            // allocate it now (the later default allocation is skipped when
            // one already exists, so this wiring is not clobbered).
            if (lmask_lev[lev].empty() || !lmask_lev[lev][0]) {
                lmask_lev[lev].resize(1);
                auto ngv = surf_mf.nGrowVect(); ngv[2] = 0;
                lmask_lev[lev][0] = std::make_unique<iMultiFab>(
                    surf_mf.boxArray(), surf_mf.DistributionMap(), 1, ngv);
            }

            // comp 0 = land-sea mask (1 = land), comp 1 = SST
            // UPSTREAM_ISSUES 29m. Copy only as many ghost cells as the
            // DESTINATION actually owns. A fresh start allocates sst_lev here
            // with surf_mf's full ghost vector (4,4,4), but the restart path
            // allocates it with ng[2]=0 (ERF_Checkpoint.cpp: "ng =
            // vars_new[lev][Vars::cons].nGrowVect(); ng[2]=0;") and this branch
            // is skipped because sst_lev[lev] is then non-empty. Copying with
            // (4,4,4) into a (4,4,0) destination addresses k=-4..4 on a fab that
            // owns one k-plane, so k<0 indexes 4*kstride = 82944 bytes BELOW the
            // fab and k>0 the same above -- straight through whatever the arena
            // put there. AMReX guards this with
            //   BL_ASSERT(dst.nGrowVect().allGE(nghost) && ...)   [MultiFab.cpp]
            // which is compiled out in Release, so it corrupts silently. It
            // killed five 24-h runs by overwriting the CACHED FillBoundary tag
            // vector of vars_new[0][cons] with SST values (~288 K), after which
            // FB_local_copy_gpu dereferenced a float as a fab pointer.
            IntVect ng_sst_copy = surf_mf.nGrowVect();
            ng_sst_copy.min(sst_lev[lev][0]->nGrowVect());
            MultiFab::Copy(*sst_lev[lev][0], surf_mf, 1, 0, 1, ng_sst_copy);

            // ---- UPSTREAM_ISSUES item 60 ----------------------------------
            // comp 0 is the LAND-SEA MASK and was allocated above but NEVER
            // POPULATED. lmask_lev stayed at its default, so LMASK_0 read a
            // single unique value of 0 -- every cell classified SEA, and
            // ERF_SurfaceLayer.cpp:845 routed every column (land included) to
            // the SST branch.
            //
            // This is the SECOND of item 60's two defects and only becomes
            // visible once the first is fixed. The first was in the reader
            // above (the unskipped lat/lon planes): while that stood, comp 0
            // was identically 0 anyway, so populating lmask from it would have
            // changed nothing -- measured, 0 land cells of 18432. Both fixes
            // are required, and neither alone is sufficient.
            //
            // lmask is an iMultiFab and surf_mf is Real, so this cannot be a
            // MultiFab::Copy -- convert per cell. Ghosts are clamped to what
            // the DESTINATION owns for exactly the reason the 29m comment
            // above gives: the restart path allocates with ng[2]=0 and an
            // over-wide copy walks off the fab silently in Release.
            {
                IntVect ng_lm = surf_mf.nGrowVect();
                ng_lm.min(lmask_lev[lev][0]->nGrowVect());
                ng_lm[2] = 0;
                for (MFIter mfi(*lmask_lev[lev][0]); mfi.isValid(); ++mfi) {
                    const Box& gbx = mfi.growntilebox(ng_lm);
                    const Array4<int>&        lm = lmask_lev[lev][0]->array(mfi);
                    const Array4<const Real>& sm = surf_mf.const_array(mfi);
                    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        lm(i,j,k) = (sm(i,j,k,0) >= Real(0.5)) ? 1 : 0;
                    });
                }
                Print() << "HindCast land-sea mask -> lmask: "
                        << lmask_lev[lev][0]->sum(0) << " land cells of "
                        << lmask_lev[lev][0]->boxArray().numPts() << std::endl;
            }

            // comp 2 = surface albedo (ERA5 fal). Register only when the
            // frames carry it (comp 2 >= 0); radiation falls back to its
            // land/sea constants otherwise. Refresh contents at every read
            // like SST (radiation keeps the raw pointer).
            if (surf_mf.max(2) >= zero) {
                if (alb_lev[lev].empty() || !alb_lev[lev][0]) {
                    alb_lev[lev].resize(1);
                    alb_lev[lev][0] = std::make_unique<MultiFab>(
                        surf_mf.boxArray(), surf_mf.DistributionMap(), 1,
                        surf_mf.nGrowVect());
                    Print() << "Hindcast surface frames provide albedo: "
                            << "per-column field registered for radiation (min/max "
                            << surf_mf.min(2) << " / " << surf_mf.max(2) << ")" << std::endl;
                }
                MultiFab::Copy(*alb_lev[lev][0], surf_mf, 2, 0, 1, surf_mf.nGrowVect());
            }

            // Sanitize: ERA5 SST carries ~9999 fill values over land, and the
            // bilinear interpolation blends them into coastal sea cells (values
            // of several thousand K were observed). qsat() of such values
            // overflows -- fatally in single precision -- and one bad surface
            // cell NaNs its whole column through the vertically-implicit solve.
            // Land cells get a benign placeholder (their surface temperature
            // comes from the land model / T0, never from SST); sea cells are
            // clamped to a physical ocean range so fill-contaminated coastal
            // cells stay finite (warm-biased in at most a 1-2 cell halo of the
            // ~25 km source data; flagged in the run documentation).
            {
                auto const& m_arrs = surf_mf.const_arrays();
                auto const& s_arrs = sst_lev[lev][0]->arrays();
                ParallelFor(*sst_lev[lev][0], sst_lev[lev][0]->nGrowVect(),
                            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    if (m_arrs[box_no](i,j,0) >= myhalf) {
                        s_arrs[box_no](i,j,k) = Real(288.0);  // land placeholder
                    } else {
                        s_arrs[box_no](i,j,k) =
                            amrex::min(amrex::max(s_arrs[box_no](i,j,k),
                                                  Real(271.0)), Real(305.0));
                    }
                });
                Gpu::streamSynchronize();
            }
            sst_lev[lev][0]->FillBoundary(geom[lev].periodicity());

            // Build TSK: the COMPLETE surface temperature field, which is what
            // SurfaceLayer reads for every cell once TSK exists (m_ignore_sst,
            // ERF_SurfaceLayer.H:465). Two sources, because the two surfaces
            // are interpolated with opposite masks:
            //   water -> sst_lev, already clamped to 271-305 by the sweep above
            //   land  -> comp 3, the land-masked skin temperature
            //
            // NOT comp 1 over land: that is the water-only interpolation, so
            // its land cells are the -1 no-valid-source flag. Copying it here
            // and clamping to a physical band turns -1 into the floor and
            // yields a COLDER constant than the soil constant it replaced --
            // measured 200 K domain-wide on the first attempt. The failure is
            // silent in every log line except the field range itself.
            //
            // A coastal land cell whose 2x2 source stencil holds no land point
            // falls back to sst_lev (288 K placeholder over land), i.e. the
            // pre-existing behaviour, rather than to a fabricated value.
            {
                IntVect ng_tsk = surf_mf.nGrowVect();
                ng_tsk.min(tsk_lev[lev][0]->nGrowVect());   // 29m: dest ghosts, not source
                ng_tsk.min(sst_lev[lev][0]->nGrowVect());

                auto const& sm_arrs = surf_mf.const_arrays();
                auto const& ss_arrs = sst_lev[lev][0]->const_arrays();
                auto const& t_arrs  = tsk_lev[lev][0]->arrays();
                ParallelFor(*tsk_lev[lev][0], ng_tsk,
                            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                {
                    const Real sea = ss_arrs[box_no](i,j,k);
                    if (sm_arrs[box_no](i,j,0) >= myhalf) {
                        const Real lnd = sm_arrs[box_no](i,j,0,3);
                        t_arrs[box_no](i,j,k) = (lnd > Real(0.0))
                            ? amrex::min(amrex::max(lnd, Real(200.0)), Real(340.0))
                            : sea;
                    } else {
                        t_arrs[box_no](i,j,k) = sea;
                    }
                });
                Gpu::streamSynchronize();
                tsk_lev[lev][0]->FillBoundary(geom[lev].periodicity());

                // Report land and sea separately: a single min/max over the
                // field cannot distinguish "land is live" from "land fell back
                // to the sea placeholder everywhere", which is the exact
                // failure this block exists to avoid.
                ReduceOps<ReduceOpMin,ReduceOpMax,ReduceOpMin,ReduceOpMax> rop;
                ReduceData<Real,Real,Real,Real> rdat(rop);
                using RT = typename decltype(rdat)::Type;
                for (MFIter mfi(*tsk_lev[lev][0]); mfi.isValid(); ++mfi) {
                    const Box& b = mfi.validbox();
                    const auto& t = tsk_lev[lev][0]->const_array(mfi);
                    const auto& m = surf_mf.const_array(mfi);
                    rop.eval(b, rdat, [=] AMREX_GPU_DEVICE (int i,int j,int k) -> RT {
                        const bool L = (m(i,j,0) >= myhalf);
                        const Real big = Real(1.e10);
                        return { L ?  t(i,j,k) :  big, L ? t(i,j,k) : -big,
                                 L ?  big : t(i,j,k),  L ? -big : t(i,j,k) };
                    });
                }
                RT hv = rdat.value(rop);
                Real lmin = get<0>(hv), lmax = get<1>(hv);
                Real wmin = get<2>(hv), wmax = get<3>(hv);
                ParallelDescriptor::ReduceRealMin(lmin);
                ParallelDescriptor::ReduceRealMax(lmax);
                ParallelDescriptor::ReduceRealMin(wmin);
                ParallelDescriptor::ReduceRealMax(wmax);
                Print() << "HindCast skin temperature -> tsk: land "
                        << lmin << " - " << lmax << " K, sea "
                        << wmin << " - " << wmax << " K" << std::endl;
            }

            auto const& mask_arrs = surf_mf.const_arrays();
            auto const& lmsk_arrs = lmask_lev[lev][0]->arrays();
            ParallelFor(*lmask_lev[lev][0], lmask_lev[lev][0]->nGrowVect(),
                        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                lmsk_arrs[box_no](i,j,k) = (mask_arrs[box_no](i,j,0) >= myhalf) ? 1 : 0;
            });
            Gpu::streamSynchronize();
            lmask_lev[lev][0]->FillBoundary(geom[lev].periodicity());
        }

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

    // Fill the time-interpolated surface state. This LinComb was commented
    // out upstream, leaving surface_state_interp UNINITIALIZED while the
    // slow-RHS/substep source terms consume it whenever hindcast_surface_bcs
    // is on (garbage land-sea mask -> unbounded bulk-coefficient tendencies,
    // NaN within one step in memory-dependent regions).
    MultiFab::LinComb(surface_state_interp[lev],
                      alpha1, surface_state_1[lev], 0,
                      alpha2, surface_state_2[lev], 0,
                      0, surface_state_interp[lev].nComp(),
                      surface_state_interp[lev].nGrow());

    // ***********************************************************************************
    // Push the ERA5 SST into the surface layer's t_surf over OCEAN points.
    //
    // Without this the sea surface is a spatially and temporally UNIFORM constant for
    // the whole run: t_surf is initialised to `default_land_surf_temp` =
    // erf.most.surf_temp (ERF_SurfaceLayer.H:445) and `set_t_surf` is only ever reached
    // from InitType::Input_Sounding (ERF.cpp:1615), so the HindCast path never touches
    // it. `fill_qsurf_with_qsat` then builds q_surf from that constant. On a domain that
    // is ~75% ocean that is a substantial deficiency in its own right, independent of
    // the tropical constants that used to sit in the bulk source term.
    //
    // Land points are left alone: over land t_surf belongs to the LSM (for this campaign
    // the MM5 soil_theta constant, since MM5 is inert), so the SST only overwrites where
    // the land-sea mask says water. The SST guard is two-sided for the same reason it is
    // in the bulk path -- this field carries fill values, and a one-sided bound lets them
    // through.
    // ***********************************************************************************
    if (m_SurfaceLayer) {
        MultiFab* tsurf = m_SurfaceLayer->get_t_surf(lev);
        if (tsurf) {
            // The SST frames carry 9999 land fills that the generator interpolates
            // against without masking, so a ring of coastal water cells comes through
            // smeared (measured on Jan-9: 10.4% of water cells fail the guard, all
            // >340 K; 0% land in the plausible band, i.e. no false accepts). Rather
            // than leave those on the deck constant -- which happened to sit 0.5 K
            // from the truth this month and would be several K out in summer, silently
            // -- fill them from valid neighbours. The clean field spans 13.1-15.7 C
            // with std 0.67 K, so neighbours are good estimates.
            // Scratch allocated ONCE and reused. Allocating these per call -- and
            // worse, allocating `nxt` inside the 12-sweep loop -- meant 12 GPU
            // MultiFab allocations every timestep, which segfaulted a 24-h run at
            // ~27.5k steps (SIGSEGV, no NaN, no assert).
            // ERF MEMBERS, not function-local statics -- see ERF.H. Statics
            // outlive amrex::Finalize() and abort at exit with CUDA 709.
            if (!m_sst_fill_sf ||
                !m_sst_fill_sf->boxArray().CellEqual(surface_state_interp[lev].boxArray())) {
                m_sst_fill_sf  = std::make_unique<MultiFab>(surface_state_interp[lev].boxArray(),
                                                            surface_state_interp[lev].DistributionMap(), 2, 1);
                m_sst_fill_nxt = std::make_unique<MultiFab>(surface_state_interp[lev].boxArray(),
                                                            surface_state_interp[lev].DistributionMap(), 2, 1);
            }
            MultiFab& sf  = *m_sst_fill_sf;
            MultiFab& nxt = *m_sst_fill_nxt;
            sf.setVal(0.0);
            for (MFIter mfi(sf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.growntilebox();
                const Array4<Real>&       a  = sf.array(mfi);
                const Array4<const Real>& ss = surface_state_interp[lev].const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real lsm = ss(i,j,0,0), sst = ss(i,j,0,1);
                    bool ok  = (lsm < Real(0.5)) && (sst > Real(271.0)) && (sst < Real(305.0));
                    a(i,j,k,0) = ok ? sst : Real(0.0);
                    a(i,j,k,1) = ok ? Real(1.0) : Real(0.0);
                });
            }
            // Iterative nearest-valid spreading: each sweep grows the valid set by one
            // cell, so 12 sweeps reach 12 cells inland from any valid water point.
            for (int it = 0; it < 12; ++it) {
                sf.FillBoundary(geom[lev].periodicity());
                MultiFab::Copy(nxt, sf, 0, 0, 2, 1);
                for (MFIter mfi(sf); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.tilebox();
                    const Array4<Real>&       b = nxt.array(mfi);
                    const Array4<const Real>& a = sf.const_array(mfi);
                    const Array4<const Real>& ss = surface_state_interp[lev].const_array(mfi);
                    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        if (a(i,j,k,1) > Real(0.5)) { return; }          // already valid
                        if (ss(i,j,0,0) >= Real(0.5)) { return; }        // land: leave to LSM
                        Real sum = Real(0.0), cnt = Real(0.0);
                        for (int dj = -1; dj <= 1; ++dj) {
                        for (int di = -1; di <= 1; ++di) {
                            if (a(i+di,j+dj,k,1) > Real(0.5)) {
                                sum += a(i+di,j+dj,k,0); cnt += Real(1.0);
                            }
                        }}
                        if (cnt > Real(0.0)) { b(i,j,k,0) = sum/cnt; b(i,j,k,1) = Real(1.0); }
                    });
                }
                MultiFab::Copy(sf, nxt, 0, 0, 2, 1);
            }

            // Push onto t_surf over water; count anything still unfilled.
            // Also score the filled cells against the originals. The fill is
            // distance-ordered (already-valid cells are skipped, so originals are
            // never modified), but a cell filled on sweep k averages neighbours that
            // were themselves filled on sweep k-1, so values flatten with distance
            // from real data. These cells are COASTAL, where the true SST gradient is
            // largest, so a collapse toward the field mean would erase land-sea
            // contrast over ~20% of the ocean.
            ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpMin,ReduceOpMax> rop2;
            ReduceData<Real,Real,Real,Real> rdat2(rop2);   // filled: sum, n, min, max
            using RT2 = typename decltype(rdat2)::Type;
            ReduceOps<ReduceOpSum,ReduceOpSum,ReduceOpMin,ReduceOpMax> rop3;
            ReduceData<Real,Real,Real,Real> rdat3(rop3);   // original: sum, n, min, max
            using RT3 = typename decltype(rdat3)::Type;
            for (MFIter mfi(*tsurf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<const Real>& a  = sf.const_array(mfi);
                const Array4<const Real>& ss = surface_state_interp[lev].const_array(mfi);
                rop2.eval(bx, rdat2, [=] AMREX_GPU_DEVICE (int i,int j,int k) -> RT2 {
                    bool wat = ss(i,j,0,0) < Real(0.5);
                    bool raw = wat && ss(i,j,0,1) > Real(271.0) && ss(i,j,0,1) < Real(305.0);
                    bool fil = wat && !raw && a(i,j,k,1) > Real(0.5);
                    Real v = fil ? a(i,j,k,0) : Real(0.0);
                    return {v, fil?Real(1.0):Real(0.0),
                            fil?a(i,j,k,0):Real(1.0e30), fil?a(i,j,k,0):Real(-1.0e30)};
                });
                rop3.eval(bx, rdat3, [=] AMREX_GPU_DEVICE (int i,int j,int k) -> RT3 {
                    bool wat = ss(i,j,0,0) < Real(0.5);
                    bool raw = wat && ss(i,j,0,1) > Real(271.0) && ss(i,j,0,1) < Real(305.0);
                    Real v = raw ? ss(i,j,0,1) : Real(0.0);
                    return {v, raw?Real(1.0):Real(0.0),
                            raw?ss(i,j,0,1):Real(1.0e30), raw?ss(i,j,0,1):Real(-1.0e30)};
                });
            }

            ReduceOps<ReduceOpSum,ReduceOpSum> reduce_op;
            ReduceData<Real,Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            for (MFIter mfi(*tsurf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.tilebox();
                const Array4<Real>&       ts = tsurf->array(mfi);
                const Array4<const Real>& a  = sf.const_array(mfi);
                const Array4<const Real>& ss = surface_state_interp[lev].const_array(mfi);
                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (ss(i,j,0,0) < Real(0.5) && a(i,j,k,1) > Real(0.5)) { ts(i,j,k) = a(i,j,k,0); }
                });
                reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    Real water = (ss(i,j,0,0) < Real(0.5)) ? Real(1.0) : Real(0.0);
                    Real raw   = (water > Real(0.0) && ss(i,j,0,1) > Real(271.0)
                                                    && ss(i,j,0,1) < Real(305.0)) ? Real(1.0) : Real(0.0);
                    Real unfilled = (water > Real(0.0) && a(i,j,k,1) < Real(0.5)) ? Real(1.0) : Real(0.0);
                    return {water - raw, unfilled};
                });
            }
            ReduceTuple hv = reduce_data.value();
            Real n_filled = amrex::get<0>(hv), n_unfilled = amrex::get<1>(hv);
            ParallelDescriptor::ReduceRealSum(n_filled);
            ParallelDescriptor::ReduceRealSum(n_unfilled);
            static bool reported = false;
            if (!reported) {
                reported = true;
                RT2 h2 = rdat2.value(); RT3 h3 = rdat3.value();
                Real fs=amrex::get<0>(h2), fn=amrex::get<1>(h2), fmin=amrex::get<2>(h2), fmax=amrex::get<3>(h2);
                Real os=amrex::get<0>(h3), on=amrex::get<1>(h3), omin=amrex::get<2>(h3), omax=amrex::get<3>(h3);
                ParallelDescriptor::ReduceRealSum(fs); ParallelDescriptor::ReduceRealSum(fn);
                ParallelDescriptor::ReduceRealMin(fmin); ParallelDescriptor::ReduceRealMax(fmax);
                ParallelDescriptor::ReduceRealSum(os); ParallelDescriptor::ReduceRealSum(on);
                ParallelDescriptor::ReduceRealMin(omin); ParallelDescriptor::ReduceRealMax(omax);
                Print() << "HindCast SST -> t_surf: " << (long)n_filled
                        << " ocean cells fell back to nearest-valid SST (frame fill/coastal smear); "
                        << (long)n_unfilled << " still unfilled (left on erf.most.surf_temp)"
                        << std::endl;
                if (fn > Real(0.0) && on > Real(0.0)) {
                    Print() << "  SST filled  cells: n=" << (long)fn << "  min " << fmin-Real(273.15)
                            << "  mean " << fs/fn-Real(273.15) << "  max " << fmax-Real(273.15) << " C\n"
                            << "  SST original cells: n=" << (long)on << "  min " << omin-Real(273.15)
                            << "  mean " << os/on-Real(273.15) << "  max " << omax-Real(273.15) << " C"
                            << std::endl;
                }
            }
        }
    }

    /*MultiFab& mf_surf_interp   = surface_state_interp[lev];

    // Fill the time-interpolated forecast states
    MultiFab::LinComb(surface_state_interp[lev],
                      alpha1, surface_state_1[lev], 0,
                      alpha2, surface_state_2[lev], 0,
                      0, mf_surf_interp.nComp(), mf_surf_interp.nGrow());

    std::string pltname = "plt_interp_surface";
    Vector<std::string> varnames_plot_mf = {"ls_mask", "SST"};

    const MultiFab& src = vars_new[0][0];

    MultiFab plot_mf(src.boxArray(),
                     src.DistributionMap(),
                     2, 0);

    plot_mf.setVal(0.0);

    for (MFIter mfi(plot_mf); mfi.isValid(); ++mfi) {
        const Array4<Real> &plot_mf_arr = plot_mf.array(mfi);
        const Array4<Real> &surf_mf_arr = surface_state_1[0].array(mfi);

        const Box& bx = mfi.validbox();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            plot_mf_arr(i,j,k,0) = surf_mf_arr(i,j,0);
            plot_mf_arr(i,j,k,1) = surf_mf_arr(i,j,1);
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
