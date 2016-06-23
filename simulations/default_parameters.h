#include "../dynein_struct.h"

double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
double T = 293.0; // K

double dt = 1e-12;

double Lt = 15.0; // nm, guess - not sure how DNA tail-bridge works
double Ls = 21.22; // nm, derived from PyMol dynein crystal struct 3VKH, 212.2 angstroms

// tail domain radius, not sure how to get since no info on DNA tail-bridge
double fake_radius_t = 1.5;  // nm
// motor domain radius, derived from PyMol, motor radius 148.6 angstroms
double fake_radius_m = 1.48; // nm
// binding domain radius, derived from PyMol, binding radius 14.78 angstroms
double fake_radius_b = 0.14; // nm

// water viscosity is about 0.7 mPa s = 7e-4 Pa s = 7e-4 m^2 * kg/s / m^3
// mu = 7e-4 kg/s / m
//  ... thus this is ...
// mu = 7e-4 kg/s / (1e9 nm)
//    = 7e-13 kg/s / nm
double water_viscosity_mu = 7e-13; // kg/(s*nm)

double gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
double gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
double gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

double ct = 0.001; // force*distance = energy = nm^2 * kg / s^2
double cm = 0.5; // ???
double cb = 0.5; // ???

double ONEBOUND_UNBINDING_FORCE = 8e12;
double BOTHBOUND_UNBINDING_FORCE = 2e12;

double MICROTUBULE_REPULSION_FORCE = 30.0; // N/nm

double MICROTUBULE_BINDING_DISTANCE = 0.2; // nm

double RAND_INIT_SEED = 0;

onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.0 * M_PI,
  0.6 * M_PI
};

bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  0.5 * M_PI,
  0.5 * M_PI,
  1.0 * M_PI,
  1.5 * M_PI,
  0.5 * M_PI
};
