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
  double bba = 0, bma = 0, fma = 0, fba = 0, bbx = 0, bby = 4; // FIXME

  // we will probably need Jin's code to output which foot steps
  // do we care about nearbound or farbound

  MTRand* rand = new MTRand(0.0); // FIXME
  Dynein_onebound* dyn_ob = new Dynein_onebound(bba, bma, fma, fba, bbx, bby,
                                                NEARBOUND,
                                                NULL,
                                                NULL,
                                                NULL,
                                                rand);

  double t = 0;
  long long iter = 0;
  bool stillStepping = true;

  while(stillStepping){

    double binding_prob = dyn_ob->get_binding_rate()*dt;
    double unbinding_prob = dyn_ob->get_unbinding_rate()*dt;

    // deal with falling off of microtubule
    if (rand->rand() < unbinding_prob) {
      delete dyn_ob;
      dyn_ob = NULL;
      printf("I fell off of the microtubule!");
      exit(1);
    }


    //
    else if (rand->rand() < binding_prob){ // bind condition
      printf("I took a step! Final L = %f", dyn_ob->get_bbx()-dyn_ob->get_ubx());
      exit(0);
    }

    else{
      t += dt; // iterate time
      iter ++;

      double old_bba = dyn_ob->get_bba();
      double old_bma = dyn_ob->get_bma();
      double old_uba = dyn_ob->get_uba();
      double old_uma = dyn_ob->get_uma();


      bool accept_step = false;
      int attempts = 0;
      const long long max_attempts = 1e6;

      while(!accept_step){
        if(attempts > max_attempts){
          printf("Over %lld attempts needed to avoid a NaN state in onebound at time %g, something must be wrong. Exiting.\n", max_attempts, t);
          exit(1);
        }
        dyn_ob->set_bba(old_bba);
        dyn_ob->set_bma(old_bma);
        dyn_ob->set_uba(old_uba);
        dyn_ob->set_uma(old_uma);

        double temp_bba = dyn_ob->get_bba() + dyn_ob->get_d_bba() * dt;
        double temp_bma = dyn_ob->get_bma() + dyn_ob->get_d_bma() * dt;
        double temp_uma = dyn_ob->get_uma() + dyn_ob->get_d_uma() * dt;
        double temp_uba = dyn_ob->get_uba() + dyn_ob->get_d_uba() * dt;

        dyn_ob->set_bba(temp_bba);
        dyn_ob->set_bma(temp_bma);
        dyn_ob->set_uma(temp_uma);
        dyn_ob->set_uba(temp_uba);

        accept_step = dyn_ob->update_velocities(); //NOTE: double check why this is a bool and not void
        attempts++;
      }
    }
  }
}
