#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "default_parameters.h"
#include "dynein_struct.h"
#include "simulation_defaults.h"


int main() {
  double bba = 0.3, bma = 0.5, fma = 0.25, fba = 0.15, bbx = 0, bby = 4; // FIXME

  // we will probably need Jin's code to output which foot steps
  // do we care about nearbound or farbound

  MTRand* rand = new MTRand(0.0); // FIXME
  Dynein_onebound* dynein = new Dynein_onebound(bba, bma, fma, fba, bbx, bby,
                                                NEARBOUND,
                                                NULL,
                                                NULL,
                                                NULL,
                                                rand);

  double t = 0;
  long long iter = 0;
  bool stillStepping = true;

  while(stillStepping){

    double binding_prob = dynein->get_binding_rate()*dt;
    if (binding_prob > 0) printf("binding_prob is %g\n", binding_prob);

    if (rand->rand() < binding_prob) {
      // We are going to bind!
      printf("I took a step! Final L = %f", dynein->get_bbx()-dynein->get_ubx());
      exit(0);
    } else {
      t += dt; // iterate time
      iter ++;

      if (binding_prob > 0) printf("binding domain at %g %g\n", dynein->get_ubx(), dynein->get_uby());

      bool accept_step = false;

      while (!accept_step) {
        double temp_bba = dynein->get_bba() + dynein->get_d_bba() * dt;
        double temp_bma = dynein->get_bma() + dynein->get_d_bma() * dt;
        double temp_uma = dynein->get_uma() + dynein->get_d_uma() * dt;
        double temp_uba = dynein->get_uba() + dynein->get_d_uba() * dt;

        dynein->set_bba(temp_bba);
        dynein->set_bma(temp_bma);
        dynein->set_uma(temp_uma);
        dynein->set_uba(temp_uba);

        accept_step = dynein->update_velocities(); //NOTE: double check why this is a bool and not void
      }
    }
  }
}
