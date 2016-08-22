#include "dynein_struct.h"

double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
double T = 293.0; // K

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

double ct = kb*310.15; // 0.5*c*<theta^2> = 0.5*kb*T
double cm = kb*310.15; // body temperature
double cb = kb*310.15;

double D = 2*kb*T / ((gt + gm + gb) / 3);

double tau = (Lt/2 + Ls/2)*(Lt/2 + Ls/2) / D;

double dt = 1e-11;

double bba_correlation_time = 1e-5; // from ob_PE_correlation_vs_time for current parameters

double binding_preexponential_factor = 1e15;
double unbinding_preexponential_factor = 1e15;

double ONEBOUND_UNBINDING_FORCE = 1e12;
double BOTHBOUND_UNBINDING_FORCE = 1e12;

double MICROTUBULE_REPULSION_FORCE = 30.0; // N/nm

double MICROTUBULE_BINDING_DISTANCE = 0.2; // nm

double RAND_INIT_SEED = 0;

onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.3 * M_PI,
  0.6 * M_PI
};

double t_nma = acos(Lt/(2*Ls));
bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  M_PI - 2*t_nma + M_PI/3,
  t_nma,
  2*t_nma + M_PI/3,
  2*M_PI - t_nma,
  2*t_nma - M_PI/3
};
