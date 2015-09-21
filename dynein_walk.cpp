#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"
int runtime;
double dt;

  /*
   * For every timestep, call update_velocities to update internal velocities.
   * Then do Euler's method to update internal coordinates and log output.
   * update_velocities must be called after every set_x command to update
   * internal velocities.
   */

void simulateProtein(Dynein* dyn, double dt, double tf) {
  double t = 0;
  
  double temp_bba;
  double temp_bma;
  double temp_fma;
  double temp_fba;

  srand(time(NULL));
  
  while( t < tf ) {
    dyn->get_unbinding_probability();
    
    if (fabs(dyn->get_fbx()) < 0.1) {
      if ((rand() % 100) / 100.0 < dyn->get_binding_probability()) {
	dyn->switch_near_far_state();
      } else if ((rand() % 100) / 100 < dyn->get_unbinding_probability()) {
	dyn->unbind();
      }
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
    
    dyn->log(t);
    
    t += dt;
  }
}


/* *********************************** MAIN ****************************************** */

int main(int argvc, char **argv) {

  if (argvc != 7) {
    printf("Error. Usage: ./walk inctime bla_init mla_init mra_init bra_init.\n");
    exit(EXIT_FAILURE);
  }
  
  dt = strtod(argv[1], NULL);
  runtime  = atoi(argv[2]);
  double bba_init = strtod(argv[3], NULL) * M_PI;
  double bma_init = strtod(argv[4], NULL) * M_PI;
  double fma_init = strtod(argv[5], NULL) * M_PI;
  double fba_init = strtod(argv[6], NULL) * M_PI;
  
  Dynein* dyn = new Dynein(bba_init, bma_init, fma_init, fba_init, // Initial angles
			   FARBOUND,                               // Initial state
			   NULL,                               // Optional custom internal forces
			   NULL,                               // Optional custom brownian forces
			   NULL);                              // Optional custom equilibrium angles
  
  dyn->resetLog();
  simulateProtein(dyn, dt, runtime);
  free(dyn);
  dyn = NULL;
  return 0;
}
