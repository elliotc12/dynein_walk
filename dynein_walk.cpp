#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "dynein_struct.h"

/*
 * For every timestep, call update_velocities to update internal velocities.
 * Then do Euler's method to update internal coordinates and log output.
 * update_velocities must be called after every set_a command to update
 * internal velocities.
 */

double dt = 1e-12; // s
double runtime;

double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
double T = 293.0; // K

double lt = 15.0;   // nm, guess - not sure how DNA tail-bridge works
double ls = 21.22; // nm, derived from PyMol dynein crystal struct 3VKH, 212.2 angstroms
double Lt = lt;
double Ls = ls; // FIXME remove ls and lt in favor of Ls and Lt

// tail domain radius, not sure how to get since no info on DNA
// tail-bridge
double fake_radius_t = 1.5;  // nm
// motor domain radius, derived from PyMol, motor radius 148.6
// angstroms
double fake_radius_m = 1.48; // nm
// binding domain radius, derived from PyMol, binding radius 14.78
// angstroms
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

/* *********************************** MAIN ****************************************** */

int main(int argvc, char **argv) {
  if (argvc != 6) {
    printf("Error. Usage: ./walk runtime bla_init mla_init mra_init bra_init.\n");
    exit(EXIT_FAILURE);
  }

  runtime = strtod(argv[1], NULL) * dt;

  MTRand* rand = new MTRand(RAND_INIT_SEED);

  //onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;

  // double bba_init = strtod(argv[2], NULL) * M_PI + eq.bba;
  // double bma_init = strtod(argv[3], NULL) * M_PI + bba_init + eq.bma - M_PI;
  // double uma_init = strtod(argv[4], NULL) * M_PI + bma_init + eq.ta;
  // double uba_init = strtod(argv[5], NULL) * M_PI + uma_init + M_PI - eq.uma;

  // Dynein_onebound* dyn_ob = new Dynein_onebound(
  // 				    bba_init, bma_init, uma_init, uba_init, // Initial angles
  // 				    0, 0,                // Starting coordinates
  // 				    FARBOUND,            // Initial state
  // 				    NULL,                // Optional custom internal forces
  // 				    NULL,                // Optional custom brownian forces
  // 				    NULL,                // Optional custom equilibrium angles
  //                                rand);
  // Dynein_bothbound *dyn_bb = 0;

  Dynein_onebound* dyn_ob = NULL;
  
  double R = sqrt(2*kb*T/(gm*dt)); // random scaling factor
  bothbound_forces no_forces = {0,0,0,0,0,0,0,0,0,0};
  bothbound_forces out_forces = {0,0,-R,0,0,0,R,0,0,0};
  
  Dynein_bothbound* dyn_bb = new Dynein_bothbound(5*M_PI/6,     // nma_init		      
						  7*M_PI/6,      // fma_init		      
						  0,             // nbx_init		      
						  0,             // nby_init		      
						  Lt,            // L -- equilateral roof  
						  &no_forces,    // internal forces	      
						  &out_forces,   // brownian forces	      
						  NULL,          // equilibrium angles     
						  rand);         // MTRand                 

  double t = 0;

  double distance_traveled = 0;
  double run_length = 0;
  double steps = 0;

  FILE* data_file = fopen("data.txt", "w+");
  FILE* run_file = fopen("run.txt", "w+");
  FILE* config_file = fopen("config.txt", "w+");

  resetLogs(data_file, config_file, runtime);
  dyn_bb->log(t, data_file);
  t += dt;

  while( t < runtime ) {
    if (dyn_ob != NULL)
      while (t < runtime) { // loop as long as it is onebound
	if (rand->rand() < dyn_ob->get_unbinding_rate()*dt) {
	  // this is the case where we fall off and become zerobound!
	  printf("unbinding.");
	  dyn_ob->log(t, data_file);
	  goto end_simulation;
	  return EXIT_SUCCESS;
	} else if (rand->rand() < dyn_ob->get_binding_rate()*dt*1e40) { // testing, bind rate huge
	  // switch to bothbound
	  steps++;
	  distance_traveled += fabs(dyn_ob->get_ubx() - dyn_ob->get_bbx());
	  run_length = (dyn_ob->get_bbx() + dyn_ob->get_ubx()) / 2.0;

	  dyn_bb = new Dynein_bothbound(dyn_ob, rand);
	  delete dyn_ob;
	  printf("switch to bothbound!\n");
	  break;
	} else { // move like normal
	  double temp_bba = dyn_ob->get_bba() + dyn_ob->get_d_bba() * dt;
	  double temp_bma = dyn_ob->get_bma() + dyn_ob->get_d_bma() * dt;
	  double temp_uma = dyn_ob->get_uma() + dyn_ob->get_d_uma() * dt;
	  double temp_uba = dyn_ob->get_uba() + dyn_ob->get_d_uba() * dt;

	  dyn_ob->set_bba(temp_bba);
	  dyn_ob->set_bma(temp_bma);
	  dyn_ob->set_uma(temp_uma);
	  dyn_ob->set_uba(temp_uba);
	}
	dyn_ob->update_velocities();
	dyn_ob->log(t, data_file);
	t += dt;
      }

    if (dyn_bb != NULL) {
      while (t < runtime) { // loop as long as it is bothbound
	bool unbind_near = rand->rand() < dyn_bb->get_near_unbinding_rate();
	bool unbind_far = rand->rand() < dyn_bb->get_far_unbinding_rate();
	if (unbind_near && unbind_far) {
	  printf("THEY BOTH WANT TO FALL OFF TOGETHER!!!\n");
	  printf("WHAT SHOULD THEY DO????\n");
	}
	if (unbind_near) {
	  // switch to farbound
	  dyn_ob = new Dynein_onebound(dyn_bb, rand, FARBOUND);
	  delete dyn_bb;
	  printf("switch to onebound!\n");
	  break;
	} else if (unbind_far) {
	  // switch to nearbound
	  dyn_ob = new Dynein_onebound(dyn_bb, rand, NEARBOUND);
	  delete dyn_bb;
	  printf("switch to onebound!\n");
	  break;
	} else {
	  // printf("d_nma: %g\n", dyn_bb->get_d_nma());
	  // printf("d_fma: %g\n", dyn_bb->get_d_fma());
	  double temp_nma = dyn_bb->get_nma() + dyn_bb->get_d_nma()*dt;
	  double temp_fma = dyn_bb->get_fma() + dyn_bb->get_d_fma()*dt;

	  dyn_bb->set_nma(temp_nma);
	  dyn_bb->set_fma(temp_fma);
	}

	dyn_bb->update_velocities();
	dyn_bb->log(t, data_file);
	t += dt;
      }
    }
  }

 end_simulation: log_run(run_file, t, run_length, distance_traveled, steps);
  fclose(data_file);
  fclose(run_file);
  fclose(config_file);

  return 0;
}
