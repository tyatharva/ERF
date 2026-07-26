#ifndef ERF_WEATHERDATAINTERPOLATION_H_
#define ERF_WEATHERDATAINTERPOLATION_H_
/**
 * Trilinear interpolation of weather forecast data onto the simulation mesh
 * The coarse weather forecast data is interpolated in time first to get the forecast
 * at the current time, and then spatially interpolated onto the simulation mesh
 */

#include <filesystem>
#include <stdexcept>
#include <fstream>
#include <cstdint>
#include <vector>
#include <limits>
#include "ERF.H"
#include "ERF_EOS.H"
#include "ERF_ReadCustomBinaryIC.H"
#include "ERF_Interpolation_Bilinear.H"

using namespace amrex;
namespace fs = std::filesystem;

enum class MultiFabType { CC, NC };

/**
 * Height (m) of the lowest frame level that carries independent data, set by
 * FillForecastStateMultiFabs and consumed by init_thermo_from_hindcast, which
 * blends from the ERA5 2-m anchor up to it. Negative until the first frame is
 * read.
 */
static amrex::Real s_frame_zlow = amrex::Real(-1.0);

/**
 * ERA5 surface anchor, loaded once and shared by the two consumers: the
 * hydrostatic initialization and (when erf.hindcast_blend_bdy_theta is set) the
 * near-surface theta blend applied to the frame itself.
 */
static amrex::Gpu::DeviceVector<amrex::Real> s_sp_d, s_zo_d, s_t2_d;
static int  s_anchor_nx = 0, s_anchor_ny = 0;
static bool s_anchor_ok = false;
static bool s_anchor_tried = false;

/**
 * Read the ERA5 surface anchor -- surface pressure, the height that pressure is
 * valid at, and the 2-m air temperature -- on the level-0 cell-centre grid.
 * Written by precip_check/make_sfc_anchor.py; see that file for the format.
 *
 * Returns false (and says why) if the file is absent, has the wrong header, or
 * is short. The caller treats that as fatal rather than falling back, because
 * the fallback is a DIFFERENT initialization and a silent revert would be
 * indistinguishable in the output from the new path having run.
 */
static bool
read_sfc_anchor (const std::string& fname, const int nx_dom, const int ny_dom,
                 Gpu::DeviceVector<Real>& sp_d,
                 Gpu::DeviceVector<Real>& zo_d,
                 Gpu::DeviceVector<Real>& t2_d)
{
    constexpr std::int32_t magic = 0x45524653;      // 'ERFS'
    const Long npts = Long(nx_dom) * Long(ny_dom);
    Vector<Real> h(3*npts, Real(0.0));
    int ok = 0;

    if (ParallelDescriptor::IOProcessor())
    {
        std::ifstream ifs(fname, std::ios::binary);
        if (!ifs.good()) {
            Print() << "HindCast IC anchor: cannot open '" << fname << "'" << std::endl;
        } else {
            std::int32_t hdr[3] = {0,0,0};
            ifs.read(reinterpret_cast<char*>(hdr), std::streamsize(3*sizeof(std::int32_t)));
            if (hdr[0] != magic || hdr[1] != nx_dom || hdr[2] != ny_dom) {
                Print() << "HindCast IC anchor: header mismatch in '" << fname
                        << "' -- got magic/nx/ny = " << hdr[0] << "/" << hdr[1] << "/" << hdr[2]
                        << ", expected " << magic << "/" << nx_dom << "/" << ny_dom << std::endl;
            } else {
                std::vector<double> buf(std::size_t(3*npts));
                const std::streamsize nbytes = std::streamsize(3*npts*sizeof(double));
                ifs.read(reinterpret_cast<char*>(buf.data()), nbytes);
                if (ifs.gcount() != nbytes) {
                    Print() << "HindCast IC anchor: short read from '" << fname << "' ("
                            << ifs.gcount() << " of " << nbytes << " bytes)" << std::endl;
                } else {
                    for (Long n = 0; n < 3*npts; ++n) { h[n] = Real(buf[std::size_t(n)]); }
                    ok = 1;
                }
            }
        }
    }
    ParallelDescriptor::Bcast(&ok, 1, ParallelDescriptor::IOProcessorNumber());
    if (!ok) { return false; }
    ParallelDescriptor::Bcast(h.data(), int(3*npts), ParallelDescriptor::IOProcessorNumber());

    sp_d.resize(std::size_t(npts)); zo_d.resize(std::size_t(npts)); t2_d.resize(std::size_t(npts));
    Gpu::copyAsync(Gpu::hostToDevice, h.begin(),        h.begin()+  npts, sp_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, h.begin()+  npts, h.begin()+2*npts, zo_d.begin());
    Gpu::copyAsync(Gpu::hostToDevice, h.begin()+2*npts, h.begin()+3*npts, t2_d.begin());
    Gpu::streamSynchronize();

    Print() << "HindCast IC anchor: read ERA5 sp / orography / t2m from '" << fname
            << "' (" << nx_dom << " x " << ny_dom << ")" << std::endl;
    return true;
}

/**
 * Load the anchor once into the file-scope statics. Returns false if no anchor
 * file was configured; aborts if one was configured but cannot be read, because
 * silently continuing would run a different initialization than the one asked for.
 */
static bool
ensure_sfc_anchor (const int nx_dom, const int ny_dom)
{
    if (s_anchor_tried) {
        return s_anchor_ok && (s_anchor_nx == nx_dom) && (s_anchor_ny == ny_dom);
    }
    s_anchor_tried = true;

    std::string fname;
    { amrex::ParmParse pp("erf"); pp.query("hindcast_sfc_anchor_file", fname); }
    if (fname.empty()) { return false; }

    s_anchor_ok = read_sfc_anchor(fname, nx_dom, ny_dom, s_sp_d, s_zo_d, s_t2_d);
    if (!s_anchor_ok) {
        Abort("erf.hindcast_sfc_anchor_file was set but could not be read. Refusing to fall "
              "back to the p_0-at-sea-level anchor silently -- the two paths produce "
              "different initial states and the difference would not be visible downstream.");
    }
    s_anchor_nx = nx_dom; s_anchor_ny = ny_dom;
    return true;
}

/**
 * erf.hindcast_blend_bdy_theta: apply the 2-m theta blend to the FRAME, so every
 * consumer of the frame sees the same near-surface profile.
 *
 * Without it the blend lives only in init_thermo_from_hindcast, so the interior is
 * initialized with the blended profile while the lateral relaxation drives the band
 * toward the raw frame -- which below the lowest frame level is the CLAMPED,
 * vertically uniform value. Measured at t = 0 on the Jan-9 case, that mismatch is
 * +1.90 K at the first cell centre in the band (d = 0..3), +1.55 K at the second,
 * decaying to +0.16 K by the fifth (z ~ 250 m) and to zero above the lowest frame
 * level. Applying it here makes the interior and the relaxation target the same
 * field by construction.
 */
/**
 * erf.hindcast_frame_from_T: rebuild the frame's theta and rho from its recovered
 * temperature and a hydrostatic pressure anchored on ERA5 sp.
 *
 * The frame's T is correct (item 21 measures it to 0.03 K at every level); its
 * theta and rho are not, both because the (rho, theta) pair it stores implies a
 * pressure ~18 hPa too low. T is exactly recoverable from that pair -- the same
 * wrong pressure cancels -- so the correct fields can be rebuilt with no new data.
 */
static bool
frame_from_T ()
{
    static const bool b = [] {
        int v = 0; amrex::ParmParse pp("erf");
        pp.query("hindcast_frame_from_T", v); return v != 0; }();
    return b;
}

static bool
blend_bdy_theta ()
{
    static const bool b = [] {
        int v = 0; amrex::ParmParse pp("erf");
        pp.query("hindcast_blend_bdy_theta", v); return v != 0; }();
    return b;
}


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

    // ------------------------------------------------------------------
    // Drop the fabricated bottom level.
    //
    // erftools emits a surface level (z = 0) whose data is a BYTE-IDENTICAL
    // copy of the level above it (z = 155.07 m for the Channel Islands case):
    // measured maxdiff = 0.000e+00 in all eight fields -- rho, uvel, vvel,
    // wvel, theta, qv, qc, qr. It is the same source level written twice under
    // two different heights, and it carries no independent information; it is
    // the same class of error as erftools' ~300 m constant downward
    // displacement of the levels that DO carry data.
    //
    // Removing it is numerically a no-op for the interpolation itself:
    // bilinear_interpolation clamps z below zvec[0], so cells under 155 m
    // received the level-1 value either way. What it changes is that nothing
    // downstream can mistake the duplicate for a surface observation. In
    // particular init_thermo_from_hindcast needs the height of the lowest
    // level that actually carries data, so it can anchor the layer beneath it
    // on ERA5 t2m rather than extend a fabricated value to the ground.
    //
    // The test is exact equality over every field, so a frame without the bug
    // (or with a genuinely distinct surface level) is left untouched.
    // ------------------------------------------------------------------
    if (zvec_h.size() >= 2 && zvec_h[0] < zvec_h[1])
    {
        const Long nxy = Long(xvec_h.size()) * Long(yvec_h.size());
        Vector<Vector<Real>*> flds = {&rho_h, &uvel_h, &vvel_h, &wvel_h,
                                      &theta_h, &qv_h, &qc_h, &qr_h};
        bool dup = true;
        for (auto* f : flds) {
            if (Long(f->size()) < 2*nxy) { dup = false; break; }
            for (Long n = 0; n < nxy; ++n) {
                if ((*f)[n] != (*f)[nxy+n]) { dup = false; break; }
            }
            if (!dup) { break; }
        }
        if (dup) {
            for (auto* f : flds) { f->erase(f->begin(), f->begin()+nxy); }
            const Real z_dropped = zvec_h[0];
            zvec_h.erase(zvec_h.begin());
            static bool announced = false;
            if (!announced) {
                announced = true;
                Print() << "HindCast frames: dropped the bottom level at z = " << z_dropped
                        << " m -- byte-identical to the level at z = " << zvec_h[0]
                        << " m in all 8 fields (erftools duplicate). Lowest level with "
                        << "independent data is now z = " << zvec_h[0] << " m." << std::endl;
            }
        }
    }
    // ------------------------------------------------------------------
    // RETRACTED: erf.hindcast_frame_z_offset (item 19f) has been removed.
    //
    // It relabelled the frame levels upward by ~305 m on the theory that erftools
    // displaced the data in height. The audit (item 21) shows it does not: the
    // frame's T, qv, u and v are correct at the labelled height, and only the
    // pressure implied by (rho, theta) is wrong. Shifting the levels made the
    // model sample air from 305 m lower -- warmer and moister -- which cancelled
    // a warm-dry bias introduced elsewhere and produced a 32x precipitation gain.
    // That gain was a real measurement of a compensating error, not a fix: with
    // the shift applied, T at 3130 m went from 272.01 K (correct, 272.04 in ERA5)
    // to 273.86 K, and qv from 0.00264 to 0.00341, i.e. +29%. It made correct
    // fields wrong.
    //
    // The pressure error is handled where it belongs, in frame_from_T() below.
    // ------------------------------------------------------------------

    s_frame_zlow = zvec_h[0];

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
            // qc and qr were being interpolated onto the ERF mesh and then
            // dropped on the floor (audit item 21, section 4). Driving with
            // hydrometeors zeroed at inflow costs 120-165 km of spin-up on this
            // domain (Roberge et al. 2024, GMD 17, 1497-1510), and the observed
            // Jan-9 event sits 144-192 km from the inflow edge -- inside it.
            if (ncomp_cons > RhoQ3_comp) {
                fine_cons_arr(i,j,k,RhoQ2_comp) = tmp_qc;
                fine_cons_arr(i,j,k,RhoQ3_comp) = tmp_qr;
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

    // ------------------------------------------------------------------
    // Rebuild the frame's theta and rho from its TEMPERATURE.
    //
    // The audit (item 21) measures the frame's T, qv, u and v as correct at the
    // labelled height -- T to within 0.03 K at every level -- while its theta is
    // +1.6 K and its rho -1.9% wrong. Both of those follow from one cause: the
    // (rho, theta) pair erftools stores implies a pressure ~18 hPa too low.
    //
    //     theta = T (p0/p)^kappa      p 1.914% low -> theta +1.57 K
    //     rho   = p /(R_d T (1+..))   p 1.914% low -> rho   -1.914%
    //
    // T is therefore EXACTLY recoverable from the pair: the same wrong p that
    // corrupted both cancels in T = p(rho,theta)/(R_d rho (1+..)), which is what
    // getTgivenRandRTh computes. No new data is needed.
    //
    // So: recover T, integrate pressure hydrostatically from ERA5's surface
    // pressure, and derive a consistent theta and rho. That DROPS the erftools
    // pressure error instead of cancelling it against another one -- which is
    // what the previous two attempts did. Item 19 paired the frame's (wrong)
    // theta with a correct surface pressure and so converted the pressure error
    // into +1.7 to +2.2 K of warming; item 19f then cancelled that warming by
    // shifting the levels 305 m, which made T and qv wrong where they had been
    // right. Both are retracted.
    //
    // Applied HERE, to the frame itself, so the initial condition, the boundary
    // planes and the relaxation target all consume the same corrected field.
    // Fixing only one of them is item 19e's defect over again.
    //
    // The near-surface anchor is folded in as a blend on T rather than on theta,
    // because t2m is a temperature; blending theta against a pressure that is
    // itself being solved for is what made the earlier version hard to reason
    // about.
    // ------------------------------------------------------------------
    if (frame_from_T() &&
        ensure_sfc_anchor(geom[lev].Domain().length(0), geom[lev].Domain().length(1)))
    {
        const Real  rdOcp   = solverChoice.rdOcp;
        const Real  l_grav  = solverChoice.gravity;
        const Real  z_low   = s_frame_zlow;
        const Real* sp_p    = s_sp_d.data();
        const Real* zo_p    = s_zo_d.data();
        const Real* t2_p    = s_t2_d.data();
        const auto  dom     = geom[lev].Domain();
        const int   nxd     = dom.length(0);
        const int   do_blend = blend_bdy_theta() ? 1 : 0;

        for (MFIter mfi(erf_mf_cons); mfi.isValid(); ++mfi) {
            const Box  gbx  = mfi.growntilebox();
            const int  klo  = gbx.smallEnd(2);
            const int  khi  = gbx.bigEnd(2);
            const Box  b2d  = makeSlab(gbx, 2, klo);
            const auto z_arr = (a_z_phys_nd) ? a_z_phys_nd->const_array(mfi) : Array4<const Real>{};
            const Array4<Real>& c_arr = erf_mf_cons.array(mfi);
            const int ncomp_c = erf_mf_cons.nComp();

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                // The anchor is defined on valid domain columns only; clamp so the
                // ghost ring uses its nearest interior column.
                const int ic = amrex::min(amrex::max(i, dom.smallEnd(0)), dom.bigEnd(0));
                const int jc = amrex::min(amrex::max(j, dom.smallEnd(1)), dom.bigEnd(1));
                const Long n2 = Long(jc)*Long(nxd) + Long(ic);
                const Real p_orog = sp_p[n2];
                const Real z_orog = zo_p[n2];
                const Real t2m    = t2_p[n2];

                const Real z_base = (z_arr(i,j,klo) + z_arr(i,j,klo+1)) / two;

                Real p_prev = p_orog;
                Real r_prev = getRhogivenTandPress(t2m, p_orog, Real(0.0));
                Real z_prev = z_orog;

                // Start at the domain floor, NOT at the box's ghost k. Below-ground
                // ghost cells have z < 0, where bilinear_interpolation returns 0 by
                // design; recovering T from rho = 0 divides by zero and the NaN then
                // propagates up the whole column through the recursion.
                for (int k = amrex::max(klo, dom.smallEnd(2)); k <= khi; ++k)
                {
                    const Real qv_k = (ncomp_c > RhoQ1_comp) ? c_arr(i,j,k,RhoQ1_comp) : Real(0.0);

                    // T recovered from the frame's own (rho, theta): the erftools
                    // pressure error cancels exactly in this combination.
                    const Real r_in = c_arr(i,j,k,Rho_comp);
                    if (!(r_in > Real(0.0))) { continue; }   // unfilled//below-ground cell
                    Real T_k = getTgivenRandRTh(r_in, r_in*c_arr(i,j,k,RhoTheta_comp), qv_k);

                    if (do_blend) {
                        const Real z_c = (z_arr(i,j,k) + z_arr(i,j,k+1)) / two;
                        Real w = (z_low > z_base) ? (z_c - z_base)/(z_low - z_base) : Real(1.0);
                        w = amrex::min(amrex::max(w, Real(0.0)), Real(1.0));
                        T_k = t2m + (T_k - t2m)*w;
                    }

                    // p(k) = p(k-1) - dz*g*(rho(k)+rho(k-1))/2 with rho = p/(R_d T_v).
                    // T is fixed here, so rho is LINEAR in p and this inverts in closed
                    // form -- no fixed point needed, unlike the theta-based version.
                    const Real dz_loc = (z_arr(i,j,k) + z_arr(i,j,k+1))/two - z_prev;
                    const Real cinv   = Real(1.0) /
                        (R_d * T_k * (Real(1.0) + (R_v/R_d)*qv_k));
                    const Real p_k = (p_prev - myhalf*dz_loc*l_grav*r_prev) /
                                     (Real(1.0) + myhalf*dz_loc*l_grav*cinv);
                    const Real r_k = p_k * cinv;

                    c_arr(i,j,k,Rho_comp)      = r_k;
                    c_arr(i,j,k,RhoTheta_comp) = getThgivenTandP(T_k, p_k, rdOcp);

                    p_prev = p_k; r_prev = r_k;
                    z_prev = (z_arr(i,j,k) + z_arr(i,j,k+1))/two;
                }
            });
        }
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
            else if (nvar==HindcastBdyVars::QC)  { src = &fcons; scomp = RhoQ2_comp;  } // plain qc
            else if (nvar==HindcastBdyVars::QR)  { src = &fcons; scomp = RhoQ3_comp;  } // plain qr
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
 * Method (erf.hindcast_sfc_anchor_file set): specify the THERMODYNAMIC profile
 * and solve for the mass field, which is the direction real.exe works in.
 * theta comes from the frame, pressure is integrated hydrostatically upward
 * from ERA5's surface pressure, and density follows from the equation of state.
 * Base state and thermodynamic state are then the same field by construction,
 * so buoyancy is exactly zero at init and the stratification is the observed
 * one. qv comes from the frame. Momenta remain zero (spin-up from rest).
 *
 * This replaces the inverse ordering -- frame density -> erf_enforce_hse ->
 * derived theta -- which discarded the frame's theta entirely and substituted
 * whatever theta the pressure integration happened to produce. Measured on
 * Jan-9 2023 that was +3.77 K at 155 m growing with height, and it inverted
 * the marine layer: theta FELL 3.35 K through the lowest 150 m over the whole
 * ocean, i.e. a domain-wide absolutely unstable surface layer at t = 0 where
 * ERA5 has a stable one. Two separate errors fed it:
 *
 *   1. erf_enforce_hse anchors p = p_0 = 101325 Pa at z = 0 in EVERY column,
 *      a fixed standard-atmosphere sea-level pressure. It cannot represent the
 *      synoptic pressure field, which is the entire content of a storm.
 *   2. The frame's own bottom level is a byte-identical duplicate of the level
 *      above it, so the lowest 155 m had no data at all -- and substituting the
 *      frame theta without fixing that gives constant rho AND constant theta
 *      below 155 m, hence constant p, hence zero layer thickness and non-finite
 *      RRTMGP optical depths at step 1.
 *
 * So the anchor is ERA5 sp (carried hydrostatically from ERA5's orography
 * height to ERF's 3-km terrain height), and the layer below the lowest frame
 * level with real data is a linear-in-z blend from the ERA5 2-m air
 * temperature up to that level. The marine layer is then stable because the
 * observations say it is, not because an integration error made it so.
 *
 * Without the knob the old path runs unchanged, gated by
 * erf.hindcast_ic_frame_theta, so every other HindCast case is untouched.
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
    if (l_has_moist) {
        // Forecast state stores PLAIN qv (velocity convention)
        MultiFab::Copy(qv_hse, fcons, RhoQ1_comp, 0, 1, ngv);
    }

    const Box& l_domain = geom[lev].Domain();
    const int nx_dom = l_domain.length(0);
    const int ny_dom = l_domain.length(1);

    const bool have_anchor = ensure_sfc_anchor(nx_dom, ny_dom);

    if (have_anchor && s_frame_zlow <= zero) {
        Abort("HindCast IC: the height of the lowest frame level is unset -- "
              "FillForecastStateMultiFabs must run before init_thermo_from_hindcast.");
    }

    if (!have_anchor)
    {
        MultiFab::Copy(r_hse, fcons, Rho_comp, 0, 1, ngv);
        erf_enforce_hse(lev, r_hse, p_hse, pi_hse, th_hse, qv_hse, z_phys_cc[lev]);
    }
    else
    {
        const Real l_gravity = solverChoice.gravity;
        const Real rdOcp     = solverChoice.rdOcp;
        const Real z_low     = s_frame_zlow;       // lowest frame level carrying data (m)
        const Real* sp_p     = s_sp_d.data();
        const Real* zo_p     = s_zo_d.data();
        const Real* t2_p     = s_t2_d.data();
        // When the blend is applied to the frame itself, f_arr already carries the
        // blended profile and applying it again here would compound it.
        const int   pre_blended = blend_bdy_theta() ? 1 : 0;
        const int   nxd      = nx_dom;
        const int   dlo_z    = l_domain.smallEnd(2);
        const int   dhi_z    = l_domain.bigEnd(2);
        const bool  lmoist   = l_has_moist;

        for (MFIter mfi(r_hse); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            if (bx.smallEnd(2) != dlo_z || bx.bigEnd(2) != dhi_z) {
                Abort("HindCast IC: the hydrostatic integration needs a whole column on one "
                      "box. Set amr.max_grid_size_z >= amr.n_cell[2].");
            }
            const Box b2d = makeSlab(bx, 2, dlo_z);

            const Array4<Real      >& r_arr  = r_hse.array(mfi);
            const Array4<Real      >& p_arr  = p_hse.array(mfi);
            const Array4<Real      >& pi_arr = pi_hse.array(mfi);
            const Array4<Real      >& th_arr = th_hse.array(mfi);
            const Array4<Real      >& qv_arr = qv_hse.array(mfi);
            const Array4<Real const>& f_arr  = fcons.const_array(mfi);
            const Array4<Real const>& zcc    = z_phys_cc[lev]->const_array(mfi);

            ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                const Long  n2     = Long(j)*Long(nxd) + Long(i);
                const Real  p_orog = sp_p[n2];      // ERA5 surface pressure  [Pa]
                const Real  z_orog = zo_p[n2];      // height it is valid at  [m]
                const Real  t2m    = t2_p[n2];      // 2-m air temperature    [K]

                // theta of the near-surface air in ERF's DRY convention; moisture
                // enters everywhere below through theta_m = theta*(1 + R_v/R_d*qv),
                // which is what getRhogivenThetaPress applies.
                const Real th_sfc = getThgivenTandP(t2m, p_orog, rdOcp);
                const Real qv_sfc = lmoist ? f_arr(i,j,dlo_z,RhoQ1_comp) : zero;

                // The theta blend runs from the FIRST CELL CENTRE, not from
                // z_orog. ERA5's 0.25 deg orography and ERF's 3-km terrain are
                // different surfaces -- near the coast ERA5 smears land elevation
                // out over water, so z_orog exceeded the first cell centre in a
                // band of ocean columns. Blending from z_orog there drove the
                // weight negative, it clamped to 0, and two or more adjacent
                // cells were all set to theta_2m: a constant-theta layer, i.e.
                // exactly the fabricated-uniform-layer defect this restructure
                // exists to remove. Measured that way: 847 ocean layers with
                // d(theta) = 0 and up to 2.4 K of departure from the frame over
                // coastal land. Anchoring on the first cell centre makes the
                // weight rise strictly monotonically with k, so d(theta) > 0
                // follows from theta_frame(z_low) > theta_2m alone.
                //
                // The PRESSURE integration still starts at (z_orog, sp), which is
                // where ERA5's surface pressure is actually valid; dz is negative
                // for the first step wherever z_orog sits above the first cell
                // centre, which correctly raises p there.
                const Real z_base = zcc(i,j,dlo_z);

                Real p_prev = p_orog;
                Real r_prev = getRhogivenThetaPress(th_sfc, p_orog, rdOcp, qv_sfc);
                Real z_prev = z_orog;

                for (int k = dlo_z; k <= dhi_z; ++k)
                {
                    const Real z_c  = zcc(i,j,k);
                    const Real qv_k = lmoist ? f_arr(i,j,k,RhoQ1_comp) : zero;

                    // Blend the 2-m anchor up to the lowest frame level that carries
                    // data. At and above z_low the weight is 1 and this is EXACTLY
                    // the frame value; below it the frame value is the clamped
                    // z_low value, so this is a linear-in-z profile running from the
                    // observed 2-m air temperature to the lowest real observation.
                    // Where the terrain already reaches above z_low the weight is 1
                    // everywhere and the frame profile is used unmodified.
                    // amrex::min/max bind by const reference, so the namespace-scope
                    // constexpr `zero`/`one` would be odr-used and are not device symbols.
                    Real w = (z_low > z_base) ? (z_c - z_base)/(z_low - z_base) : Real(1.0);
                    w = amrex::min(amrex::max(w, Real(0.0)), Real(1.0));
                    if (pre_blended) { w = Real(1.0); }
                    const Real th_k = th_sfc + (f_arr(i,j,k,RhoTheta_comp) - th_sfc)*w;

                    // p(k) = p(k-1) - dz*g*(rho(k) + rho(k-1))/2 with rho = rho(p,theta):
                    // the same trapezoidal rule erf_enforce_hse uses, so the base state
                    // sits in the model's OWN discrete hydrostatic balance rather than
                    // merely a continuous one. rho depends on p only as p^(1/Gamma), so
                    // the fixed point contracts by dz*g*rho/(Gamma*p) ~ 4e-3 per sweep;
                    // 4 sweeps is exact to single precision.
                    //
                    // p is strictly decreasing with height by construction here
                    // (dz > 0, g > 0, rho > 0), which is what the layer-thickness
                    // requirement in the radiation driver actually needs.
                    const Real dz_loc = z_c - z_prev;
                    Real p_k = p_prev - dz_loc*l_gravity*r_prev;
                    Real r_k = r_prev;
                    for (int it = 0; it < 4; ++it) {
                        r_k = getRhogivenThetaPress(th_k, p_k, rdOcp, qv_k);
                        p_k = p_prev - dz_loc*l_gravity*myhalf*(r_k + r_prev);
                    }
                    r_k = getRhogivenThetaPress(th_k, p_k, rdOcp, qv_k);

                    r_arr (i,j,k) = r_k;
                    p_arr (i,j,k) = p_k;
                    pi_arr(i,j,k) = getExnergivenP(p_k, rdOcp);
                    th_arr(i,j,k) = th_k;
                    qv_arr(i,j,k) = qv_k;

                    p_prev = p_k; r_prev = r_k; z_prev = z_c;
                }

                // klo-1 ghost: hydrostatic extension downward at constant theta,
                // matching what erf_enforce_hse leaves there. physbcs_base runs
                // straight afterwards and resets it; filled for parity.
                {
                    const int  km     = dlo_z - 1;
                    const Real dz_loc = zcc(i,j,dlo_z) - zcc(i,j,km);
                    const Real th_k   = th_arr(i,j,dlo_z);
                    const Real qv_k   = qv_arr(i,j,dlo_z);
                    const Real p_0z   = p_arr(i,j,dlo_z);
                    const Real r_0z   = r_arr(i,j,dlo_z);
                    Real p_k = p_0z + dz_loc*l_gravity*r_0z;
                    for (int it = 0; it < 4; ++it) {
                        const Real rr = getRhogivenThetaPress(th_k, p_k, rdOcp, qv_k);
                        p_k = p_0z + dz_loc*l_gravity*myhalf*(rr + r_0z);
                    }
                    r_arr (i,j,km) = getRhogivenThetaPress(th_k, p_k, rdOcp, qv_k);
                    p_arr (i,j,km) = p_k;
                    pi_arr(i,j,km) = getExnergivenP(p_k, rdOcp);
                    th_arr(i,j,km) = th_k;
                    qv_arr(i,j,km) = qv_k;
                }
            });
        }
         r_hse.FillBoundary(geom[lev].periodicity());
         p_hse.FillBoundary(geom[lev].periodicity());
        pi_hse.FillBoundary(geom[lev].periodicity());
        th_hse.FillBoundary(geom[lev].periodicity());
        qv_hse.FillBoundary(geom[lev].periodicity());
    }

    (*physbcs_base[lev])(base_state[lev], 0, base_state[lev].nComp(), base_state[lev].nGrowVect());

    static const int l_ic_frame_theta = [] {
        int v=0; amrex::ParmParse pp("erf");
        pp.query("hindcast_ic_frame_theta", v); return v; }();
    // With the anchor path, th_hse ALREADY holds the frame theta (blended to the
    // 2-m anchor near the ground), so reading th_arr below is reading the frame.
    // The gate only chooses between the two OLD behaviours.
    const int ic_frame_theta = have_anchor ? 0 : l_ic_frame_theta;

    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const Box& gbx = mfi.growntilebox(1);
        const Array4<Real      >& cons_arr = cons.array(mfi);
        const Array4<Real const>& r_arr    = r_hse.const_array(mfi);
        const Array4<Real const>& th_arr   = th_hse.const_array(mfi);
        const Array4<Real const>& f_arr    = fcons.const_array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            cons_arr(i,j,k,Rho_comp)      = r_arr(i,j,k);
            // theta comes from the INTERPOLATED FRAME, not from the hydrostatic base
            // state. This previously read th_arr (= th_hse, the theta erf_enforce_hse
            // derives while integrating dp/dz from the interpolated rho), so the frame's
            // theta was read, horizontally interpolated, vertically interpolated onto the
            // ERF levels and then DISCARDED -- only rho and qv survived into the
            // thermodynamic state, and the base state's hydrostatic error became the
            // state's error.
            //
            // Measured on Jan-9 2023: model theta 293.07 K at 155 m against the frame's
            // 289.30 K at the same height, +3.77 K, growing with height (equivalent to a
            // 721 m displacement at 150 m rising to 824 m at 589 m) -- the signature of
            // integration error accumulating upward. It also produced a superadiabatic
            // lowest 150 m over the whole ocean (theta falling 3.35 K with height) where
            // ERA5 has a stable marine layer, i.e. a domain-wide absolutely unstable
            // surface layer at t = 0.
            //
            // Note the lateral forcing was NEVER affected: the boundary planes take theta
            // from fcons (line ~523), so the boundaries relaxed toward the correct frame
            // theta while the interior started ~4 K warmer.
            // GATED (erf.hindcast_ic_frame_theta, default 0 = old behaviour).
            // Substituting the frame theta alone is NOT sufficient: the frame's bottom
            // level duplicates level 1, so rho and theta are both constant below 155 m
            // and p = p0(rho R theta/p0)^gamma comes out CONSTANT over the lowest five
            // ERF levels -> zero layer thickness -> RRTMGP produces non-finite optical
            // depths and the run aborts at step 1. The old th_hse path masked this
            // because erf_enforce_hse made theta vary hydrostatically, at the cost of
            // being +3.77 K wrong. The correct fix is to integrate hydrostatically FROM
            // the frame theta (specify the thermodynamic profile, solve for the mass
            // field), which is a restructure, not a substitution.
            cons_arr(i,j,k,RhoTheta_comp) = (ic_frame_theta)
                                          ? r_arr(i,j,k) * f_arr(i,j,k,RhoTheta_comp)
                                          : r_arr(i,j,k) * th_arr(i,j,k);
            if (l_has_moist) {
                cons_arr(i,j,k,RhoQ1_comp) = r_arr(i,j,k) * f_arr(i,j,k,RhoQ1_comp);
                // Seed cloud and rain water too, so the interior does not have to
                // grow condensate the frame already knows about.
                if (cons_arr.nComp() > RhoQ3_comp && f_arr.nComp() > RhoQ3_comp) {
                    cons_arr(i,j,k,RhoQ2_comp) = r_arr(i,j,k) * f_arr(i,j,k,RhoQ2_comp);
                    cons_arr(i,j,k,RhoQ3_comp) = r_arr(i,j,k) * f_arr(i,j,k,RhoQ3_comp);
                }
            }
        });
    }

    Print() << "HindCast init: base state and thermodynamic state rebuilt from "
            << "the interpolated ERA5 frame (theta/qv coupling); lev " << lev
            << " rho min/max " << cons.min(Rho_comp) << " " << cons.max(Rho_comp) << std::endl;

    if (have_anchor)
    {
        // The two structural properties this restructure exists to guarantee are
        // cheap to check here, so check them in the log rather than only in the
        // plotfile: p strictly decreasing with height (equivalently a non-zero
        // pressure thickness for every layer, which is what the radiation driver
        // needs), and a stable surface layer.
        const int dlo_z = l_domain.smallEnd(2);
        const int dhi_z = l_domain.bigEnd(2);

        ReduceOps<ReduceOpMin,ReduceOpMax,ReduceOpMin,ReduceOpMin> reduce_op;
        ReduceData<Real,Real,Real,Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(p_hse); mfi.isValid(); ++mfi) {
            const Box  bx  = mfi.validbox();
            const Box  b2d = makeSlab(bx, 2, dlo_z);
            const Array4<Real const>& p  = p_hse.const_array(mfi);
            const Array4<Real const>& th = th_hse.const_array(mfi);
            reduce_op.eval(b2d, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int) -> ReduceTuple
            {
                Real dp_min = std::numeric_limits<Real>::max();
                Real dth_min = std::numeric_limits<Real>::max();
                for (int k = dlo_z; k < dhi_z; ++k) {
                    dp_min  = amrex::min(dp_min , p (i,j,k) - p (i,j,k+1));
                    dth_min = amrex::min(dth_min, th(i,j,k+1) - th(i,j,k));
                }
                return { p(i,j,dlo_z), p(i,j,dlo_z), dp_min, dth_min };
            });
        }
        ReduceTuple hv = reduce_data.value(reduce_op);
        Real psfc_min = amrex::get<0>(hv), psfc_max = amrex::get<1>(hv);
        Real dp_min   = amrex::get<2>(hv), dth_min  = amrex::get<3>(hv);
        ParallelDescriptor::ReduceRealMin(psfc_min);
        ParallelDescriptor::ReduceRealMax(psfc_max);
        ParallelDescriptor::ReduceRealMin(dp_min);
        ParallelDescriptor::ReduceRealMin(dth_min);

        Print() << "HindCast init [anchor]: lowest frame level with data z = " << s_frame_zlow
                << " m\n    p(k=0) min/max = " << psfc_min << " / " << psfc_max << " Pa"
                << "\n    min layer pressure thickness = " << dp_min
                << " Pa  (must be > 0)"
                << "\n    min d(theta) across any layer = " << dth_min
                << " K  (negative = superadiabatic somewhere; expected over heated land, "
                << "not over the ocean)" << std::endl;
        if (dp_min <= zero) {
            Abort("HindCast IC: non-positive layer pressure thickness after hydrostatic "
                  "integration -- the radiation driver cannot run on this state.");
        }
    }

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
