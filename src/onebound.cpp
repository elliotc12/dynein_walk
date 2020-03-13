#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>

#include <iostream>

#include "default_parameters.h"
#include "simulation_defaults.h"
#include "dynein_struct.h"

// see dynein_struct.h for possible command line parameters (extern double ... )
// they are all defined in default_parameters.h

void get_command_line_flags(int argc, char** argv, double *bba, double *bma, double *uma, double *uba, double *state, double *k, double *sticky_rate){
  // for(int i=1; i<argc; i++){
  //   fprintf(stderr, "%s ", argv[i]);
  // }
  // fprintf(stderr, "\n");
  assert(argc==25);
  int i = 1;
  low_affinity_binding_rate = atof(argv[i++]);
  *sticky_rate = atof(argv[i++]);
  cb = atof(argv[i++]);
  cm = atof(argv[i++]);
  ct = atof(argv[i++]);
  Ls = atof(argv[i++]);
  Lt = atof(argv[i++]);
  // r_t
  fake_radius_t = atof(argv[i++]);
  // r_m
  fake_radius_m = atof(argv[i++]);
  // r_b
  fake_radius_b = atof(argv[i++]);
  // seed
  RAND_INIT_SEED = atof(argv[i++]);
  dt = atof(argv[i++]);
  // eqb
  onebound_post_powerstroke_internal_angles.bba = atof(argv[i]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.nba = atof(argv[i])* M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.fba = atof(argv[i++]) * M_PI / 180.0;
  // eqmpre
  onebound_post_powerstroke_internal_angles.uma = atof(argv[i++]) * M_PI / 180.0;
  // eqmpost
  onebound_post_powerstroke_internal_angles.bma = atof(argv[i]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.nma = atof(argv[i]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.fma = atof(argv[i++]) * M_PI / 180.0;
  // eqt
  onebound_post_powerstroke_internal_angles.ta = atof(argv[i]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.ta = atof(argv[i++]) * M_PI / 180.0;
  //force
  tail_force = atof(argv[i++]) * 0.6022 / atp_in_kJ_per_mol; // conversion for our force units: 1 (dG ATP kJ / mol / nm) = atp_in_kJ_per_mol * 1e-11 / 6.022 N
  // exp_unbinding_const;
  exponential_unbinding_angle_constant = atof(argv[i++]);
  *bba = atof(argv[i++]);
  *bma = atof(argv[i++]);
  *uma = atof(argv[i++]);
  *uba = atof(argv[i++]);
  *state = atof(argv[i++]);
  *k = atof(argv[i++]);
  assert(i==25);
  // note we are not using runtime as onebound runs indefinitely until binding event
}

int main(int argc, char** argv) {
  double bba, bma, uma, uba, bbx = 0, bby = 0, state, k, sticky_rate;
  get_command_line_flags(argc, argv, &bba, &bma, &uma, &uba, &state, &k, &sticky_rate);
  State s = NEARBOUND;
  // fprintf(stderr, "angles: %g %g %g %g\n state: %g\n", bba, bma, uma, uba, state);
  if (state == 1) {
    s = FARBOUND;
  }
  binding_mode = EXPONENTIAL_UNBINDING;

  Rand* rand = new Rand(RAND_INIT_SEED);

  // Does this onebound resemble with the bothbound orientation

  // fprintf(stderr, "angles: %g %g %g %g\n state: %i\n", bba, bma, uma, uba, s);
  Dynein_onebound* dynein = new Dynein_onebound(bba, bma, uma, uba, bbx, bby,
                                                s,
                                                NULL,
                                                NULL,
                                                NULL,
                                                rand);
  // fprintf(stderr, "Starting with %g %g %g %g\n",
  //         dynein->get_bba(), dynein->get_bma(), dynein->get_uma(), dynein->get_uba());

  double t = 0;
  long iter = 0;
  bool stillStepping = true;
  double cumulative_prob = 0;
  // fprintf(stderr, "I am initially %g %g\n", dynein->get_bbx(), dynein->get_bby());
  // fprintf(stderr, "I am initially %g %g\n", dynein->get_ubx(), dynein->get_uby());

  // FILE *log_f = NULL;
  // if (k == 1) {
  //   log_f = fopen("raw-data-1.dat", "w");
  // } else if (k == 3) {
  //   log_f = fopen("raw-data-3.dat", "w");
  // }
  const double sticky_prob = sticky_rate*dt;
  bool am_sticky = false;
  while(stillStepping){

    // We are sticky if we were already sticky or if we become sticky.
    am_sticky = am_sticky || rand->rand() < sticky_prob;

    double binding_prob = dynein->get_binding_rate()*dt;
    cumulative_prob += binding_prob;

    // fprintf(stderr, "sticky rate: %g  sticky prob: %g\n", sticky_rate, sticky_prob);


    //fprintf(stderr, "The time is %g\n", t);
    // if (binding_prob > 0) fprintf(stderr, "binding_prob is %g at time %g with angle %g (total %g)\n",
    //                               binding_prob, t, dynein->get_uba(), cumulative_prob);

    t += dt;  // iterate time
    iter ++;

    double old_bba = dynein->get_bba();
    double old_bma = dynein->get_bma();
    double old_uma = dynein->get_uma();
    double old_uba = dynein->get_uba();
    // fprintf(stderr, "energy: %g at time %g with cumulative %g\n", dynein->get_PE(), t, cumulative_prob);


    bool accept_step = false;
    int attempts = 0;
    while(!accept_step){
      if (attempts > 0) {
        dynein->set_bba(old_bba);
        dynein->set_bma(old_bma);
        dynein->set_uma(old_uma);
        dynein->set_uba(old_uba);
        dynein->update_velocities();    // This update_velocities() was added (?)
      }

      double temp_bba = dynein->get_bba() + dynein->get_d_bba() * dt;
      double temp_bma = dynein->get_bma() + dynein->get_d_bma() * dt;
      double temp_uma = dynein->get_uma() + dynein->get_d_uma() * dt;
      double temp_uba = dynein->get_uba() + dynein->get_d_uba() * dt;

      dynein->set_bba(temp_bba);
      dynein->set_bma(temp_bma);
      dynein->set_uma(temp_uma);
      dynein->set_uba(temp_uba);

      accept_step = dynein->update_velocities(); //NOTE: double check why this is a bool and not void

      attempts++;
    }

    if (am_sticky && binding_prob > 0 && rand->rand() < binding_prob) {
      // We are going to bind!
      // fprintf(stderr, "I took a step after %ld! Final L = %f\n =====> %.15g %.15g %.15g %.15g\n",
      //         iter, dynein->get_ubx()-dynein->get_bbx(),
      //         dynein->get_bmy(), dynein->get_ty(), dynein->get_umy(), dynein->get_uby());
      printf("{\n  'L': %g,\n  't': %g,\n}\n",
              dynein->get_ubx()-dynein->get_bbx(), t);
      // printf("L: %g,\nt: %g\n", dynein->get_bbx()-dynein->get_ubx(), t); // YAML version
      exit(0);
    }
  }
}
