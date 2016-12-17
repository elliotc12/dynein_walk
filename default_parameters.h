#include "dynein_struct.h"

// double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
const double atp_in_kJ_per_mol = 30.5; // energy stored in 1 ATP in kJ/mol
const double eV_per_kJ_per_mol = .01036410; // 1 kJ/mol = .01 eV
const double kb_eV = 8.61733034e-5; // eV/K
double kb = atp_in_kJ_per_mol*eV_per_kJ_per_mol*kb_eV; // kB in ATP energies per K
double T = 310.15; // K

double Lt = 7.0;  // nm, guess - not sure how DNA tail-bridge works
double Ls = 12.0; // nm, derived from PyMol dynein crystal struct 3VKH, 212.2 angstroms

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

double d_tail_theta = 1.5;
double d_motor_theta = 0.75;
double d_binding_theta = 0.5;

const double binding_energy_high_affinity_kJ_mol = 71; // kJ/mol
double binding_energy_high_affinity_atp = binding_energy_high_affinity_kJ_mol / atp_in_kJ_per_mol;

/* double ct = kb*310.15/d_tail_theta/d_tail_theta;   // 0.5*c*<theta^2> = 0.5*kb*T */
/* double cm = kb*310.15/d_motor_theta/d_motor_theta; // body temperature */
/* double cb = kb*310.15/d_binding_theta/d_binding_theta; */

double ct =  1*binding_energy_high_affinity_atp; // ct = 0.5 cb
double cm = 6*binding_energy_high_affinity_atp; // cm = 3*cb
double cb = 2*binding_energy_high_affinity_atp; // see thesis_stuff.pdf 'Estimating cb spring constant'

double D = 2*kb*T / ((gt + gm + gb) / 3);

double tau = (Lt/2 + Ls/2)*(Lt/2 + Ls/2) / D;

double dt = 1e-11;

double low_affinity_binding_rate = 180; //s^-1
double low_affinity_unbinding_rate = 5000; //s^-1
//double high_affinity_unbinding_rate = ; //s^-1
//double high_affinity_binding_rate = ; //s^-1

double e = exp(1.0);

double DELTA_G_FORMATION_BINDING = 1e-10;

double ONEBOUND_UNBINDING_FORCE = 1e12;
double BOTHBOUND_UNBINDING_FORCE = 1e12;

double MICROTUBULE_REPULSION_FORCE = 30.0; // N/nm
double MICROTUBULE_BINDING_DISTANCE = 0.2; // nm

double RAND_INIT_SEED = 0;

onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
  126.0 * M_PI / 180.0,
  164.0 * M_PI / 180.0,
  0.0   * M_PI / 180.0,
  74.0 * M_PI / 180.0
};

bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  126.0 * M_PI / 180.0,
  164.0 * M_PI / 180.0,
  0.0   * M_PI / 180.0,
  164.0 * M_PI / 180.0,
  126.0 * M_PI / 180.0
};

/* double t_nma = acos(Lt/(2*Ls)); */
/* bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = { */
/*   M_PI - 2*t_nma + M_PI/3, */
/*   t_nma, */
/*   2*t_nma + M_PI/3, */
/*   2*M_PI - t_nma, */
/*   2*t_nma - M_PI/3 */
/* }; */
