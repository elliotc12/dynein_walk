#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"
  /*
   * For every timestep, call update_velocities to update internal velocities.
   * Then do Euler's method to update internal coordinates and log output.
   * update_velocities must be called after every set_x command to update
   * internal velocities.
   */

extern const double dt;
double runtime;

void simulateProtein(Dynein* dyn, double tf) {
  double t = 0;
  
  double temp_bba;
  double temp_bma;
  double temp_fma;
  double temp_fba;

  srand(time(NULL));
  FILE* data_file = fopen("data.txt", "a+");
  
  while( t < tf ) {
    if ((rand() % 100) / 100 < dyn->get_unbinding_probability()) {
        dyn->unbind();
	dyn->log(t, data_file);
	exit(EXIT_SUCCESS);
    } else if ((rand() % 100) / 100 < dyn->get_binding_probability()) {
        dyn->switch_to_bothbound();
    }
    
    dyn->update_velocities();
    
    temp_bba = dyn->get_bba() + dyn->get_d_bba() * dt;
    temp_bma = dyn->get_bma() + dyn->get_d_bma() * dt;
    temp_fma = dyn->get_fma() + dyn->get_d_fma() * dt;
    temp_fba = dyn->get_fba() + dyn->get_d_fba() * dt;
    
    dyn->set_bba(temp_bba);
    dyn->set_bma(temp_bma);
    dyn->set_fma(temp_fma);
    dyn->set_fba(temp_fba);
    
    dyn->log(t, data_file);
    
    t += dt;
  }
  
  fclose(data_file);
}


/* *********************************** MAIN ****************************************** */

int main(int argvc, char **argv) {

  if (argvc != 6) {
    printf("Error. Usage: ./walk runtime bla_init mla_init mra_init bra_init.\n");
    exit(EXIT_FAILURE);
  }

  runtime  = strtod(argv[1], NULL) * dt;

  equilibrium_angles eq = near_farbound_post_powerstroke_internal_angles;

  double bba_init = strtod(argv[2], NULL) * M_PI + eq.bba;
  double bma_init = strtod(argv[3], NULL) * M_PI + bba_init + eq.ba - M_PI;
  double fma_init = strtod(argv[4], NULL) * M_PI + bma_init + eq.ta;
  double fba_init = strtod(argv[5], NULL) * M_PI + fma_init + M_PI - eq.fa;
  
  Dynein* dyn = new Dynein(bba_init, bma_init, fma_init, fba_init, // Initial angles
			   FARBOUND,                               // Initial state
			   NULL,                               // Optional custom internal forces
			   NULL,                               // Optional custom brownian forces
			   NULL);                              // Optional custom equilibrium angles
  
  dyn->resetLog();
  simulateProtein(dyn, runtime);
  free(dyn);
  dyn = NULL;
  return 0;
}
