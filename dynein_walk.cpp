#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"

/*
 * For every timestep, call update_velocities to update internal velocities.
 * Then do Euler's method to update internal coordinates and log output.
 * update_velocities must be called after every set_a command to update
 * internal velocities.
 */

extern const double dt;
double runtime;

void simulateProtein(void* dyn, double tf) {
  double t = 0;

  FILE* data_file = fopen("data.txt", "a+");
  FILE* run_file = fopen("run.txt", "a+");

  MTRand* rand = new MTRand(RAND_INIT_SEED);

  while( t < tf ) {
             // onebound case
    if (((Dynein_onebound*) dyn)->get_state() != BOTHBOUND) { // fix type
             // unbind
      if (rand->rand() < ((Dynein_onebound*) dyn)->get_unbinding_rate()*dt) {
	printf("unbinding.");
	((Dynein_onebound*) dyn)->log_run(t, run_file);
	((Dynein_onebound*) dyn)->log(t, data_file);
	exit(EXIT_SUCCESS);
      }
      else if (rand->rand() < ((Dynein_onebound*) dyn)->get_binding_rate()*dt) {
	     // switch to bothbound
	Dynein_bothbound* new_dynein
	  = new Dynein_bothbound(((Dynein_onebound*) dyn), rand);
	free(((Dynein_onebound*) dyn));
	dyn = new_dynein;
      }
      else { // move like normal
	double temp_bba,temp_bma, temp_uma, temp_uba;
	temp_bba = ((Dynein_onebound*) dyn)->get_bba() + ((Dynein_onebound*) dyn)->get_d_bba() * dt;
	temp_bma = ((Dynein_onebound*) dyn)->get_bma() + ((Dynein_onebound*) dyn)->get_d_bma() * dt;
	temp_uma = ((Dynein_onebound*) dyn)->get_uma() + ((Dynein_onebound*) dyn)->get_d_uma() * dt;
	temp_uba = ((Dynein_onebound*) dyn)->get_uba() + ((Dynein_onebound*) dyn)->get_d_uba() * dt;

	((Dynein_onebound*) dyn)->set_bba(temp_bba);
	((Dynein_onebound*) dyn)->set_bma(temp_bma);
	((Dynein_onebound*) dyn)->set_uma(temp_uma);
	((Dynein_onebound*) dyn)->set_uba(temp_uba);
      }

      ((Dynein_onebound*) dyn)->update_velocities();
      ((Dynein_onebound*) dyn)->log(t, data_file);
    }
    else {
             // bothbound case
      if (rand->rand() < ((Dynein_bothbound*) dyn)->get_near_unbinding_rate()*dt) {
	     // switch to farbound
	Dynein_onebound* new_dynein
	  = new Dynein_onebound(((Dynein_bothbound*) dyn), rand, FARBOUND);
	free(((Dynein_bothbound*) dyn));
	dyn = new_dynein;
      }
      else if (rand->rand() < ((Dynein_bothbound*) dyn)->get_far_unbinding_rate()*dt) {
	     // switch to nearbound
	Dynein_onebound* new_dynein
	  = new Dynein_onebound(((Dynein_bothbound*) dyn), rand, NEARBOUND);
	free(((Dynein_bothbound*) dyn));
	dyn = new_dynein;
      }
      else {
	double temp_nma, temp_fma;
	temp_nma = ((Dynein_bothbound*) dyn)->get_nma() + ((Dynein_bothbound*) dyn)->get_d_nma()*dt;
	temp_fma = ((Dynein_bothbound*) dyn)->get_fma() + ((Dynein_bothbound*) dyn)->get_d_fma()*dt;

	((Dynein_bothbound*) dyn)->set_nma(temp_nma);
	((Dynein_bothbound*) dyn)->set_fma(temp_fma);
      }

      ((Dynein_bothbound*) dyn)->update_velocities();
      ((Dynein_bothbound*) dyn)->log(t, data_file);
    }

    t += dt;
  }

  ((Dynein_onebound*) dyn)->log_run(tf, run_file); // fix type

  fclose(data_file);
}


/* *********************************** MAIN ****************************************** */

int main(int argvc, char **argv) {
  if (argvc != 6) {
    printf("Error. Usage: ./walk runtime bla_init mla_init mra_init bra_init.\n");
    exit(EXIT_FAILURE);
  }

  runtime  = strtod(argv[1], NULL) * dt;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;

  double bba_init = strtod(argv[2], NULL) * M_PI + eq.bba;
  double bma_init = strtod(argv[3], NULL) * M_PI + bba_init + eq.bma - M_PI;
  double uma_init = strtod(argv[4], NULL) * M_PI + bma_init + eq.ta;
  double uba_init = strtod(argv[5], NULL) * M_PI + uma_init + M_PI - eq.uma;

  Dynein_onebound* dyn = new Dynein_onebound(
				    bba_init, bma_init, uma_init, uba_init, // Initial angles
				    0, 0,                // Starting coordinates
				    FARBOUND,            // Initial state
				    NULL,                // Optional custom internal forces
				    NULL,                // Optional custom brownian forces
				    NULL);               // Optional custom equilibrium angles

  resetLog();
  simulateProtein(dyn, runtime);
  delete dyn;
  dyn = NULL;
  return 0;
}
