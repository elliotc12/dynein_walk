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

#include <iostream>

#include "default_parameters.h"
#include "simulation_defaults.h"
#include "dynein_struct.h"

// see dynein_struct.h for possible command line parameters (extern double ... )
// they are all defined in default_parameters.h


void get_command_line_flags(int argc, char** argv);

int main(int argc, char** argv) {
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
  get_command_line_flags(argc, argv);

  int temp;
  std::cin>>temp;


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

      double old_bba = dynein->get_bba();
      double old_bma = dynein->get_bma();
      double old_uba = dynein->get_uba();
      double old_uma = dynein->get_uma();
      if (binding_prob > 0 && rand->rand() > binding_prob) {
        printf("binding domain at %g %g\n", dynein->get_ubx(), dynein->get_uby());
        printf("all done\n");
        exit(0);
      }


      bool accept_step = false;
      while(!accept_step){
        dynein->set_bba(old_bba);
        dynein->set_bma(old_bma);
        dynein->set_uba(old_uba);
        dynein->set_uma(old_uma);

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

void get_command_line_flags(int argc, char** argv){
  for(int i=1; i<argc; i++){
    printf("%s ", argv[i]);
  }
  assert(argc==19);
  low_affinity_binding_rate = atof(argv[1]);
  low_affinity_unbinding_rate = atof(argv[2]);
  cb = atof(argv[3]);
  cm = atof(argv[4]);
  ct = atof(argv[5]);
  Ls = atof(argv[6]);
  Lt = atof(argv[7]);
  // r_t
  fake_radius_t = atof(argv[8]);
  // r_m
  fake_radius_m = atof(argv[9]);
  // r_b
  fake_radius_b = atof(argv[10]);
  // seed
  RAND_INIT_SEED = atof(argv[11]);
  dt = atof(argv[12]);
  // eqb
  onebound_post_powerstroke_internal_angles.bba = atof(argv[13]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.nba = atof(argv[13])* M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.fba = atof(argv[13]) * M_PI / 180.0;
  // eqmpre
  onebound_post_powerstroke_internal_angles.uma = atof(argv[14]) * M_PI / 180.0;
  // eqmpost
  onebound_post_powerstroke_internal_angles.bma = atof(argv[15]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.nma = atof(argv[15]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.fma = atof(argv[15]) * M_PI / 180.0;
  // eqt
  onebound_post_powerstroke_internal_angles.ta = atof(argv[16]) * M_PI / 180.0;
  bothbound_pre_powerstroke_internal_angles.ta = atof(argv[16]) * M_PI / 180.0;
  //force
  tail_force = atof(argv[17]) * 0.6022 / atp_in_kJ_per_mol; // conversion for our force units: 1 (dG ATP kJ / mol / nm) = atp_in_kJ_per_mol * 1e-11 / 6.022 N
  // exp_unbinding_const;
  exponential_unbinding_angle_constant = atof(argv[18]);

  // note we are not using runtime as onebound runs indefinitely until binding event
}















