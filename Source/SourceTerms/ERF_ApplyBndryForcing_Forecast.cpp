#include <AMReX_MultiFab.H>
#include <ERF_SrcHeaders.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
ApplyBndryForcing_Forecast (
  const SolverChoice& solverChoice,
  const Geometry geom,
  const Box& tbx,
  const Box& tby,
  const Box& tbz,
  const Array4<const Real>& z_phys_nd,
  const Array4<Real>& rho_u_rhs,
  const Array4<Real>& rho_v_rhs,
  const Array4<Real>& rho_w_rhs,
  const Array4<const Real>& rho_u,
  const Array4<const Real>& rho_v,
  const Array4<const Real>& rho_w,
  const Array4<const Real>& rho_u_initial_state,
  const Array4<const Real>& rho_v_initial_state,
  const Array4<const Real>& rho_w_initial_state,
  const Array4<const Real>& cons_initial_state)

{
    // Domain cell size and real bounds
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();

    amrex::ignore_unused(rho_w_initial_state);

    // Domain valid box
    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0) + 1;
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1) + 1;

    Real hindcast_lateral_sponge_length   = solverChoice.hindcast_lateral_sponge_length;
    Real hindcast_lateral_sponge_strength = solverChoice.hindcast_lateral_sponge_strength;

    bool hindcast_zhi_sponge_damping      = solverChoice.hindcast_zhi_sponge_damping;
    Real hindcast_zhi_sponge_length       = solverChoice.hindcast_zhi_sponge_length;
    Real hindcast_zhi_sponge_strength     = solverChoice.hindcast_zhi_sponge_strength;

    Real xlo_sponge_end   = ProbLoArr[0] + hindcast_lateral_sponge_length;
    Real xhi_sponge_start = ProbHiArr[0] - hindcast_lateral_sponge_length;
    Real ylo_sponge_end   = ProbLoArr[1] + hindcast_lateral_sponge_length;
    Real yhi_sponge_start = ProbHiArr[1] - hindcast_lateral_sponge_length;
    Real zhi_sponge_start = ProbHiArr[2] - hindcast_zhi_sponge_length;

    AMREX_ALWAYS_ASSERT(xlo_sponge_end   > ProbLoArr[0]);
    AMREX_ALWAYS_ASSERT(xhi_sponge_start < ProbHiArr[0]);
    AMREX_ALWAYS_ASSERT(ylo_sponge_end   > ProbLoArr[1]);
    AMREX_ALWAYS_ASSERT(yhi_sponge_start < ProbHiArr[1]);

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);

        Real x = ProbLoArr[0] + ii * dx[0];
        Real y = ProbLoArr[1] + (jj+myhalf) * dx[1];

        Real rho_u_sponge = rho_u_initial_state(i,j,k)*cons_initial_state(i,j,k,0);
        // Use the MAX band weight rather than summing the four band terms:
        // additive stacking double-forces the band-overlap corners, which
        // was measured driving standing convergence there (density growing
        // ~0.3 kg/m3 per hour at all four corners until EOS blowup ~4.6 h).
        Real xi = zero;
        if (x < xlo_sponge_end)   { xi = amrex::max(xi, (xlo_sponge_end - x)   / hindcast_lateral_sponge_length); }
        if (x > xhi_sponge_start) { xi = amrex::max(xi, (x - xhi_sponge_start) / hindcast_lateral_sponge_length); }
        if (y < ylo_sponge_end)   { xi = amrex::max(xi, (ylo_sponge_end - y)   / hindcast_lateral_sponge_length); }
        if (y > yhi_sponge_start) { xi = amrex::max(xi, (y - yhi_sponge_start) / hindcast_lateral_sponge_length); }
        if (xi > zero) {
            rho_u_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_u(i, j, k) - rho_u_sponge);
        }
    });


    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);

        Real x = ProbLoArr[0] + (ii+myhalf) * dx[0];
        Real y = ProbLoArr[1] + jj * dx[1];

        Real rho_v_sponge    = rho_v_initial_state(i,j,k)*cons_initial_state(i,j,k,0);

        // MAX band weight (see x-momentum kernel above for the rationale)
        Real xi = zero;
        if (x < xlo_sponge_end)   { xi = amrex::max(xi, (xlo_sponge_end - x)   / hindcast_lateral_sponge_length); }
        if (x > xhi_sponge_start) { xi = amrex::max(xi, (x - xhi_sponge_start) / hindcast_lateral_sponge_length); }
        if (y < ylo_sponge_end)   { xi = amrex::max(xi, (ylo_sponge_end - y)   / hindcast_lateral_sponge_length); }
        if (y > yhi_sponge_start) { xi = amrex::max(xi, (y - yhi_sponge_start) / hindcast_lateral_sponge_length); }
        if (xi > zero) {
            rho_v_rhs(i, j, k) -= hindcast_lateral_sponge_strength * xi * xi * (rho_v(i, j, k) - rho_v_sponge);
        }
    });

    ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        Real z = z_phys_nd(i,j,k);

        if(hindcast_zhi_sponge_damping){
            if (z > zhi_sponge_start) {
                Real xi = (z - zhi_sponge_start) / hindcast_zhi_sponge_length;
                rho_w_rhs(i, j, k) -= hindcast_zhi_sponge_strength * xi * xi * (rho_w(i, j, k) - zero);
            }
        }
    });
}

// Lateral sponge for the thermodynamic state: relax (rho theta) and, with a
// moisture model active, (rho qv) toward the time-interpolated forecast in
// the same quadratic-weight bands the momentum sponge uses. The forecast
// state stores PLAIN theta / qv (velocity convention); rho-weighted targets
// are formed here. Top sponge intentionally stays w-only (standard practice:
// damp reflections, don't force thermodynamics at the lid).
void
ApplyBndryForcingCC_Forecast (
  const SolverChoice& solverChoice,
  const Geometry geom,
  const Box& bx,
  const Array4<Real>& cell_rhs,
  const Array4<const Real>& cell_data,
  const Array4<const Real>& cons_initial_state,
  const bool has_moisture)
{
    auto dx = geom.CellSizeArray();
    auto ProbHiArr = geom.ProbHiArray();
    auto ProbLoArr = geom.ProbLoArray();

    const Box& domain = geom.Domain();
    int domlo_x = domain.smallEnd(0);
    int domhi_x = domain.bigEnd(0);
    int domlo_y = domain.smallEnd(1);
    int domhi_y = domain.bigEnd(1);

    Real hindcast_lateral_sponge_length   = solverChoice.hindcast_lateral_sponge_length;
    Real hindcast_lateral_sponge_strength = solverChoice.hindcast_lateral_sponge_strength;

    Real xlo_sponge_end   = ProbLoArr[0] + hindcast_lateral_sponge_length;
    Real xhi_sponge_start = ProbHiArr[0] - hindcast_lateral_sponge_length;
    Real ylo_sponge_end   = ProbLoArr[1] + hindcast_lateral_sponge_length;
    Real yhi_sponge_start = ProbHiArr[1] - hindcast_lateral_sponge_length;

    const bool relax_qv = has_moisture && (cell_data.nComp() > RhoQ1_comp) &&
                          (cons_initial_state.nComp() > RhoQ1_comp);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
        int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);

        Real x = ProbLoArr[0] + (ii+myhalf) * dx[0];
        Real y = ProbLoArr[1] + (jj+myhalf) * dx[1];

        Real xi = zero;
        if (x < xlo_sponge_end)   { xi = amrex::max(xi, (xlo_sponge_end - x)   / hindcast_lateral_sponge_length); }
        if (x > xhi_sponge_start) { xi = amrex::max(xi, (x - xhi_sponge_start) / hindcast_lateral_sponge_length); }
        if (y < ylo_sponge_end)   { xi = amrex::max(xi, (ylo_sponge_end - y)   / hindcast_lateral_sponge_length); }
        if (y > yhi_sponge_start) { xi = amrex::max(xi, (y - yhi_sponge_start) / hindcast_lateral_sponge_length); }

        if (xi > zero) {
            Real coeff = hindcast_lateral_sponge_strength * xi * xi;
            // Weight the forecast theta/qv targets by the MODEL density
            // (WRF-nudging convention: nudge the uncoupled variable).
            // Two alternatives were measured and BOTH destabilize:
            //  - direct mass-equation relaxation toward rho_f: NaN in ~60
            //    steps (acoustic substepping does not tolerate a slow mass
            //    source);
            //  - forecast-rho weighting (rho_f*theta_f targets, momentum-
            //    sponge convention): NaN in ~50 steps.
            // Known limitation of the model-rho form: nothing anchors mass,
            // and the overlap corners of the sponge bands hold standing
            // convergence (rho reached 1.92 kg/m3 at the S-E corner surface
            // cell after ~4.6 h with the unbounded-TKE binary). Mitigate at
            // the deck level if it recurs (numerical diffusion / corner
            // geometry), not with a mass source here.
            Real rho_m = cell_data(i,j,k,Rho_comp);
            Real rt_f  = rho_m * cons_initial_state(i,j,k,RhoTheta_comp);
            cell_rhs(i,j,k,RhoTheta_comp) -= coeff * (cell_data(i,j,k,RhoTheta_comp) - rt_f);
            if (relax_qv) {
                Real rq_f = rho_m * cons_initial_state(i,j,k,RhoQ1_comp);
                cell_rhs(i,j,k,RhoQ1_comp) -= coeff * (cell_data(i,j,k,RhoQ1_comp) - rq_f);
            }
        }
    });
}
