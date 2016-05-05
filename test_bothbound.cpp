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

extern const double dt;
double runtime;

/* *********************************** MAIN ****************************************** */

int main(int argvc, char **argv) {
  MTRand* rand = new MTRand(RAND_INIT_SEED);

  // typedef struct { double nba, nma, ta, fma, fba; } bothbound_equilibrium_angles;
  {
    bothbound_equilibrium_angles eq = {
      M_PI/2, M_PI/3, 0, M_PI/3, M_PI/2
    };
    Dynein_bothbound dyn_bb(M_PI/2,
                            M_PI/3,
                            0,
                            M_PI/3,
                            0,
                            NULL, NULL, &eq, rand);
    dyn_bb.update_velocities();
    printf("d nma/dt = %g\n", dyn_bb.get_d_nma());
    printf("d fma/dt = %g\n", dyn_bb.get_d_fma());
  }

  return 0;
}
