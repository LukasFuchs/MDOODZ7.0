#include "math.h"
#include "mdoodz.h"

// ============================================================================
// ThermalConvection2D_Ref — CitcomCU reference model (Fuchs & Becker 2019,
// no strain-dependent weakening), as ThermalConvection.jl models/03_IH_ref
//
// Box 8D x D (D = 2890 km), Ra = 1e5 (with eta0 at T* = 1), purely internally
// heated, Q = 10 (nondim.), no shear/adiabatic heating.
//
// Rheology: Arrhenius-like viscosity (Tackley 2000, cstv = 2)
//   eta = eta0 exp( E/(T* + 1) - E/2 ),  E = 23.03,  T* = (T - T_s)/DeltaT >= 0
// plus constant yield stress sigma_y = 4e3 (nondim.) -> C = 48.65 MPa, phi = 0,
// with a viscoplastic floor eta_vp = 1e-2 eta0 (as eta_vp in run_IH_ref.jl).
//
// BCs:  free slip on all walls
//       T = T_s on top (Dirichlet), zero heat flux at the bottom and the sides
// Initial T: T* = 1 - z* + A sin(pi z*) cos(pi x*), z* = height above the
//       bottom / D, x* = x / D (4 wavelengths over the box), A = user2
//
// user0 = T_s [C] (surface temperature), user1 = DeltaT [K], user2 = A
// ============================================================================

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  return 0;
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  const double T_s    = (input->model.user0 + zeroC) / input->scaling.T;
  const double DeltaT = input->model.user1 / input->scaling.T;
  const double A      = input->model.user2;
  const double D      = input->model.zmax - input->model.zmin;
  const double z_nd   = (coordinates.z - input->model.zmin) / D;   // 0 at the bottom, 1 at the top
  const double x_nd   = (coordinates.x - input->model.xmin) / D;

  double T_nd = 1.0 - z_nd + A * sin(M_PI * z_nd) * cos(M_PI * x_nd);
  if (T_nd < 0.0) T_nd = 0.0;

  return T_s + DeltaT * T_nd;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->materials.rho[phase];
}

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E) {
    bc.value = 0.0;
    bc.type  = 0;      // no penetration
  } else if (position == S || position == N ||
             position == SW || position == SE ||
             position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;     // free slip: dVx/dz = 0
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == S || position == N) {
    bc.value = 0.0;
    bc.type  = 0;      // no penetration
  } else if (position == W || position == E ||
             position == SW || position == SE ||
             position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;     // free slip: dVz/dx = 0
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

char SetBCPType(MdoodzInput *input, POSITION position) {
  return position == NW ? 0 : -1;   // pressure pinned in the top-left cell (P ~ 0 at the surface)
}

SetBC SetBCT(MdoodzInput *input, POSITION position, Coordinates coordinates, double gridTemperature) {
  SetBC bc;
  const double T_s = (input->model.user0 + zeroC) / input->scaling.T;
  if (position == N || position == NE || position == NW) {
    bc.value = T_s;
    bc.type  = 1;      // Dirichlet, surface temperature
  } else {
    bc.value = 0.0;
    bc.type  = 0;      // zero heat flux: bottom and sides
  }
  return bc;
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetTemperature = SetTemperature,
                  .SetDensity     = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx,
                  .SetBCVz    = SetBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
          },
  };
  RunMDOODZ("ThermalConvection2D_Ref.txt", &setup);
}
