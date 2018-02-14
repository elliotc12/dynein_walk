#include "dynein_struct.h"

// double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
const double atp_in_kJ_per_mol = 30.5; // energy stored in 1 ATP in kJ/mol
const double eV_per_kJ_per_mol = .01036410; // 1 kJ/mol = .01 eV
const double kb_eV = 8.61733034e-5; // eV/K
double kb = kb_eV / (atp_in_kJ_per_mol * eV_per_kJ_per_mol); // kB in ATP energies per K
double T = 310.15; // K

/* double Lt = 11.15; // nm, guess - not sure how DNA tail-bridge works */
double Lt = 10; // nm, guess - not sure how DNA tail-bridge works
double Ls = 22.1; // nm, derived from PyMol dynein crystal struct 3VKH, 212.2 angstroms

// tail domain radius, derived from PyMol, see thesis_stuff
double fake_radius_t = 2.16;  // nm
// motor domain radius, derived from PyMol, see thesis_stuff
double fake_radius_m = 7.36; // nm
// binding domain radius, derived from PyMol, see thesis_stuff
double fake_radius_b = 1.57; // nm

// water viscosity is about 0.7 mPa s = 7e-4 Pa s = 7e-4 m^2 * kg/s / m^3
// mu = 7e-4 kg/s / m
//  ... thus this is ...
// mu = 7e-4 kg/s / (1e9 nm)
//    = 7e-13 kg/s / nm
double water_viscosity_mu = 7e-13; // kg/(s*nm)

double gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
double gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
double gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

const double binding_energy_high_affinity_kJ_mol = 71; // kJ/mol
double binding_energy_high_affinity_atp = binding_energy_high_affinity_kJ_mol / atp_in_kJ_per_mol;

double ct = 0.1*binding_energy_high_affinity_atp; // ct = 0.5 cb
double cm = 0.1*binding_energy_high_affinity_atp; // cm = 3*cb
double cb = 0.1*binding_energy_high_affinity_atp; // see thesis_stuff.pdf 'Estimating cb spring constant', then experimental tweaking to get the 0.1s

double dt = 1e-10;

double low_affinity_binding_rate = 180; //s^-1
double low_affinity_unbinding_rate = 460; //s^-1

double exponential_unbinding_angle_constant = 1; // radian^-1
//double high_affinity_unbinding_rate = ; //s^-1
//double high_affinity_binding_rate = ; //s^-1

/* double binding_fraction = 1e-6; */

double MICROTUBULE_REPULSION_FORCE = 0.0; // N/nm
double MICROTUBULE_BINDING_DISTANCE = 0.01; // nm

//double REBINDING_IMMUNITY_TIME = 1e-8; // s
double REBINDING_IMMUNITY_TIME = 0; // s

double RAND_INIT_SEED = 0;

TRANSITION_MODES binding_mode = EXPONENTIAL_UNBINDING;

bool am_only_writing_on_crash = false;
double stepping_movie_framerate = 1e-10;

/* onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = { */
/*    63.5 * M_PI / 180.0, */
/*   136.0 * M_PI / 180.0, */
/*     0.0 * M_PI / 180.0, */
/*   160.0 * M_PI / 180.0 */
/* }; */

/* bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = { */
/*    63.5 * M_PI / 180.0, */
/*   136.0 * M_PI / 180.0, */
/*     0.0 * M_PI / 180.0, */
/*   136.0 * M_PI / 180.0, */
/*    63.5 * M_PI / 180.0 */
/* }; */

onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
   116.5 * M_PI / 180.0,
  224.0 * M_PI / 180.0,
    0.0 * M_PI / 180.0,
  200.0 * M_PI / 180.0
};

bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  116.5 * M_PI / 180.0,
  224.0 * M_PI / 180.0,
    0.0 * M_PI / 180.0,
  224.0 * M_PI / 180.0,
  116.5 * M_PI / 180.0
};
