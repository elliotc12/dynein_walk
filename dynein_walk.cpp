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

void simulateProtein(Dynein_onebound* dyn1, double tf) {
  double t = 0;
  Dynein_bothbound *dynB = 0;

  FILE* data_file = fopen("data.txt", "a+");
  FILE* run_file = fopen("run.txt", "a+");

  MTRand* rand = new MTRand(RAND_INIT_SEED);

  while( t < tf ) {
    while (t < tf) { // loop as long as it is onebound
      if (rand->rand() < dyn1->get_unbinding_rate()*dt) {
        // this is the case where we fall off and become zerobound!
        printf("unbinding.");
        dyn1->log_run(t, run_file);
        dyn1->log(t, data_file);
        return;
      } else if (rand->rand() < dyn1->get_binding_rate()*dt) {
        // switch to bothbound
        dynB = new Dynein_bothbound(dyn1, rand);
        delete dyn1;
        break;
      } else { // move like normal
        double temp_bba = dyn1->get_bba() + dyn1->get_d_bba() * dt;
        double temp_bma = dyn1->get_bma() + dyn1->get_d_bma() * dt;
        double temp_uma = dyn1->get_uma() + dyn1->get_d_uma() * dt;
        double temp_uba = dyn1->get_uba() + dyn1->get_d_uba() * dt;

        dyn1->set_bba(temp_bba);
        dyn1->set_bma(temp_bma);
        dyn1->set_uma(temp_uma);
        dyn1->set_uba(temp_uba);
      }
      dyn1->update_velocities();
      dyn1->log(t, data_file);
      t += dt;
    }
    while (t < tf) { // loop as long as it is bothbound
      bool unbind_near = rand->rand() < dynB->get_near_unbinding_rate()*dt;
      bool unbind_far = rand->rand() < dynB->get_far_unbinding_rate()*dt;
      if (unbind_near && unbind_far) {
        printf("THEY BOTH WANT TO FALL OFF TOGETHER!!!\n");
        printf("WHAT SHOULD THEY DO????\n");
      }
      if (unbind_near) {
        // switch to farbound
	dyn1 = new Dynein_onebound(dynB, rand, FARBOUND);
	delete dynB;
        break;
      } else if (unbind_far) {
        // switch to nearbound
	dyn1 = new Dynein_onebound(dynB, rand, NEARBOUND);
	delete dynB;
        break;
      } else {
	double temp_nma = dynB->get_nma() + dynB->get_d_nma()*dt;
	double temp_fma = dynB->get_fma() + dynB->get_d_fma()*dt;

	dynB->set_nma(temp_nma);
	dynB->set_fma(temp_fma);
      }

      dynB->update_velocities();
      dynB->log(t, data_file);
      t += dt;
    }

  }

  // dyn1->log_run(tf, run_file); // fix type

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
  return 0;
}
