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
  
  double temp_bla;
  double temp_mla;
  double temp_mra;
  double temp_bra;

  srand(time(NULL));
  
  while( t < tf ) {

    if (std::abs(dyn->get_bry()) < 0.1 && rand() % 2 == 0) {
      double temp_brx = dyn->get_brx();
      double temp_bry = dyn->get_bry();
      double temp_bla = dyn->get_bla();
      double temp_mla = dyn->get_mla();
      double temp_mra = dyn->get_mra();
      double temp_bra = dyn->get_bra();
      dyn->set_blx(temp_brx);
      dyn->set_bly(temp_bry);
      dyn->set_bla(temp_bra);
      dyn->set_mla(temp_mra);
      dyn->set_mra(temp_mla);
      dyn->set_bra(temp_bla);
      if (dyn->get_state() == LEFTBOUND) dyn->set_state(RIGHTBOUND);
      else dyn->set_state(LEFTBOUND);
    }
    
    dyn->update_velocities();
    
    temp_bla = dyn->get_bla() + dyn->get_d_bla() * dt;
    temp_mla = dyn->get_mla() + dyn->get_d_mla() * dt;
    temp_mra = dyn->get_mra() + dyn->get_d_mra() * dt;
    temp_bra = dyn->get_bra() + dyn->get_d_bra() * dt;
    
    dyn->set_bla(temp_bla);
    dyn->set_mla(temp_mla);
    dyn->set_mra(temp_mra);
    dyn->set_bra(temp_bra);
    
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

  //forces no_forces =    {0,0,0,0,0,0,0,0,0,0}; // blx, bly, mlx, mly, ...
  
  dt = strtod(argv[1], NULL);
  runtime  = atoi(argv[2]);
  double bla_init = strtod(argv[3], NULL) * M_PI;
  double mla_init = strtod(argv[4], NULL) * M_PI;
  double mra_init = strtod(argv[5], NULL) * M_PI;
  double bra_init = strtod(argv[6], NULL) * M_PI;
  
  Dynein* dyn = new Dynein(bla_init, mla_init, mra_init, bra_init, // Initial angles
			   RIGHTBOUND,                             // Initial state
			   NULL,                                   // Optional custom internal forces
			   NULL,                                   // Optional custom brownian forces
			   NULL);                                  // Optional custom equilibrium angles
  
  dyn->resetLog();
  simulateProtein(dyn, dt, runtime);
  free(dyn);
  dyn = NULL;
  return 0;
}
