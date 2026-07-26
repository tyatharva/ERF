#include <AMReX_MultiFab.H>
#include <ERF_SrcHeaders.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
ApplySurfaceTreatment_BulkCoeff_Mom (
  const Box& tbx,
  const Box& tby,
  const Array4<Real>& rho_u_rhs,
  const Array4<Real>& rho_v_rhs,
  const Array4<const Real>& rho_u,
  const Array4<const Real>& rho_v,
  const Array4<const Real>& cons_state,
  const Array4<const Real>& z_phys_nd,
  const Array4<const Real>& surface_state_arr)
{
    int ndrag = 1;
    Real Cd_sea = Real(0.001);
    Real Cd_land = Real(0.01);

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        if(k <= ndrag) {
            Real dz = z_phys_nd(i,j,1) - z_phys_nd(i,j,0);
            Real ls_mask = (surface_state_arr(i-1,j,0)+surface_state_arr(i,j,0))/two;
            Real weight = one - Real(k)/ndrag;  // linear decrease
            Real fac = k==0? one : zero;
            Real Cd = fac*Cd_sea*(one-ls_mask) + Cd_land*ls_mask;
            Real rho_for_u = (cons_state(i-1,j,k,0)+cons_state(i,j,k,0))/two;
            Real rho_for_v = (cons_state(i,j-1,k,0)+cons_state(i,j,k,0))/two;
            Real uvel = rho_u(i,j,k)/rho_for_u;
            Real vvel = rho_v(i,j,k)/rho_for_v;
            Real velmag = std::sqrt(uvel*uvel + vvel*vvel);
            rho_u_rhs(i, j, k) += -one*weight*Cd*velmag*rho_for_u*uvel/dz;
        }
    });


    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
       if(k <= ndrag) {
            Real dz = z_phys_nd(i,j,1) - z_phys_nd(i,j,0);
            Real ls_mask = (surface_state_arr(i,j-1,0)+surface_state_arr(i,j,0))/two;
            Real fac = k==0? one : zero;
            Real Cd = fac*Cd_sea*(one-ls_mask) + Cd_land*ls_mask;
            Real weight = one - Real(k)/ndrag;  // linear decrease
            Real rho_for_u = (cons_state(i-1,j,k,0)+cons_state(i,j,k,0))/two;
            Real rho_for_v = (cons_state(i,j-1,k,0)+cons_state(i,j,k,0))/two;
            Real uvel = rho_u(i,j,k)/rho_for_u;
            Real vvel = rho_v(i,j,k)/rho_for_v;
            Real velmag = std::sqrt(uvel*uvel + vvel*vvel);
            rho_v_rhs(i, j, k) += -one*weight*Cd*velmag*rho_for_v*vvel/dz;
        }
    });
}

void
ApplySurfaceTreatment_BulkCoeff_CC (const Box& bx,
                         const Array4<Real>& cell_rhs,
                         const Array4<const Real>& cons_state,
                         const Array4<const Real>& u_arr,
                         const Array4<const Real>& v_arr,
                         const Array4<const Real>& z_phys_cc,
                         const Array4<const Real>& surface_state_arr)
{
     ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        if(k == 0) {
            Real dz = z_phys_cc(i,j,1)-z_phys_cc(i,j,0);
            Real ls_mask = surface_state_arr(i,j,0);
            Real Ch = Real(0.0015)*(one-ls_mask);
            Real Ce = Real(0.0015)*(one-ls_mask);

            Real rho = cons_state(i,j,k, Rho_comp);
            Real rhotheta = cons_state(i,j,k, RhoTheta_comp);
            Real rhoqv = cons_state(i,j,k, RhoQ1_comp);
            Real theta = rhotheta/rho;
            Real qv = rhoqv/rho;
            Real temp = getTgivenRandRTh(rho, rhotheta, qv);
            // Cell-centered speed from the face velocities. The previous code
            // read cons components 1/2 (RhoTheta, RhoKE) as momenta, giving
            // velmag ~ theta ~ 300 "m/s" and bulk fluxes ~100x too strong --
            // strong enough to destabilize the marine boundary layer within
            // minutes of model time (verified by bisect: runs are stable with
            // this treatment disabled and blew up in the ocean block with it).
            Real uvel = myhalf*(u_arr(i,j,k) + u_arr(i+1,j,k));
            Real vvel = myhalf*(v_arr(i,j,k) + v_arr(i,j+1,k));
            Real velmag = std::sqrt(uvel*uvel + vvel*vvel);
            // Drive toward the ACTUAL sea surface, not hardcoded tropical constants.
            //
            // This previously used dT = max(0, 301 K - temp) and dq = max(0, 0.024 - qv),
            // i.e. 28 C and a saturation mixing ratio for ~28 C, while consulting
            // surface_state only for the land/sea mask -- the ERA5 SST sitting in
            // component 1 was never used.  Measured on Jan-9 2023 (a SoCal winter
            // storm): ocean k=0 had T = 294.7 K and qv = 0.01698 kg/kg against a
            // q_sat of 0.01621, i.e. the near-surface air was ALREADY supersaturated
            // with respect to its own temperature, yet the term was still injecting
            // 9.2 mm/day of moisture.  The companion 301 K heat term added 6.3 K of
            // spurious heating, which raised q_sat and let the column hold that
            // moisture until it rained out.
            //
            // The max(0,...) clamps are also gone.  They made the flux one-way -- the
            // sea could only moisten and warm the air, never dry or cool it.  Against
            // a hardcoded tropical target that was a safety net; against a real SST it
            // is a rectifier that biases moisture upward, and condensation onto a
            // cooler sea surface plus downward sensible heat flux are both physical
            // and both matter in a winter storm.
            Real sst = surface_state_arr(i,j,0,1);
            Real dT  = zero, dq = zero;
            // Guard BOTH sides.  A lower bound alone is not enough: a large fill
            // value passes it and then esat = 611.2*exp(17.67*sstC/(sstC+243.5))
            // saturates the exponential at ~4.7e7, giving qsat ~ 1e10 and an
            // instant NaN.  Measured: with a one-sided guard both the NSCBC and
            // the Davies run went NaN in every conserved component at step 1.
            if (sst > Real(200.0) && sst < Real(340.0)) {
                Real psfc  = getPgivenRTh(rhotheta, qv);
                Real sstC  = sst - Real(273.15);
                Real esat  = Real(611.2) * std::exp(Real(17.67)*sstC / (sstC + Real(243.5)));
                Real qsat  = Real(0.622) * esat / max(psfc - esat, Real(1.0));
                dT = sst  - temp;
                dq = qsat - qv;
            }
            cell_rhs(i, j, k, RhoTheta_comp) += (theta/(Real(1005.0)*temp))*rho*Ch*velmag*dT/dz;
            cell_rhs(i, j, k, RhoQ1_comp) += rho*Ce*velmag*dq/dz;
        }
    });
}
