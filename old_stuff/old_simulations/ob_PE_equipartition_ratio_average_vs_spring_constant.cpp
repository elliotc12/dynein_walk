#include <cassert>
#include <stdio.h>
#include <string.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main(int argc, char** argv) {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  int iterations = 1e5;
  T = 1000;

  Lt = 15;
  Ls = 15;

  fake_radius_t = 1.5;
  fake_radius_m = 1.5;
  fake_radius_b = 1.5;
  gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
  gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
  gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

  const char* bba_eq_title = "Bound binding";
  const char* bma_eq_title = "Bound motor";
  const char* ta_eq_title =  "Tail domain";
  const char* uma_eq_title = "Unbound motor";

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }
  char* f_appended_name = argv[1];
  char bba_eq_fname[200];
  char bma_eq_fname[200];
  char ta_eq_fname[200];
  char uma_eq_fname[200];
  char config_eq_fname[200];

  strcpy(bba_eq_fname, "data/bba_pe_equipartition_ratio_vs_c_");
  strcpy(bma_eq_fname, "data/bma_pe_equipartition_ratio_vs_c_");
  strcpy(ta_eq_fname, "data/ta_pe_equipartition_ratio_vs_c_");
  strcpy(uma_eq_fname, "data/uma_pe_equipartition_ratio_vs_c_");
  strcpy(config_eq_fname, "data/config_pe_equipartition_ratio_vs_c_");

  strcat(bba_eq_fname, f_appended_name);
  strcat(bma_eq_fname, f_appended_name);
  strcat(ta_eq_fname, f_appended_name);
  strcat(uma_eq_fname, f_appended_name);
  strcat(config_eq_fname, f_appended_name);

  strcat(bba_eq_fname, ".txt");
  strcat(bma_eq_fname, ".txt");
  strcat(ta_eq_fname, ".txt");
  strcat(uma_eq_fname, ".txt");
  strcat(config_eq_fname, ".txt");

  prepare_data_file(bba_eq_title, bba_eq_fname);
  prepare_data_file(bma_eq_title, bma_eq_fname);
  prepare_data_file(ta_eq_title,  ta_eq_fname);
  prepare_data_file(uma_eq_title, uma_eq_fname);

  const int seeds[] = {0, 1, 2};
  int seed_len = sizeof(seeds) / sizeof(int);

  // int spring_len = 100;
  // double spring_constants[spring_len];
  // double min_spring = 0.1*kb*T;
  // double max_spring = 5;
  // double d_spring = (max_spring - min_spring) / spring_len;
  // double spring = min_spring;
  // int i = 0;

  // while(spring < max_spring) {
  //   spring_constants[i] = spring;
  //   spring += d_spring;
  //   i++;
  // }

  double spring_constants[] = {0.1*kb*T, 0.3*kb*T, 0.7*kb*T,
			       kb*T, 3*kb*T, 7*kb*T,
			       10*kb*T, 30*kb*T, 70*kb*T,
			       100*kb*T, 300*kb*T, 700*kb*T,
			       1000*kb*T};
  int spring_len = sizeof(spring_constants) / sizeof(double);

  char run_msg[512];
  const char* run_msg_base = "springconst calc (";

  onebound_data eq_data;
  double eq_bb, eq_bm, eq_t, eq_um;
  eq_data.len = 1;
  eq_data.bb = &eq_bb;
  eq_data.bm = &eq_bm;
  eq_data.t =  &eq_t;
  eq_data.um = &eq_um;

  onebound_forces unused_forces;
  generic_data unused_data;
  unused_data.data = &unused_forces;

  for (int i = 0; i < spring_len; i++) {
    ct = spring_constants[i];
    cm = spring_constants[i];
    cb = spring_constants[i];

    char buf[50];
    strcpy(run_msg, run_msg_base);
    sprintf(buf, "c = %.1g, temp = %.1g, ", spring_constants[i], T);
    strcat(run_msg, buf);

    get_onebound_equipartition_ratio
      (&eq_data, &unused_data, iterations, seeds, seed_len, run_msg);

    append_data_to_file(&spring_constants[i], &eq_bb, 1, bba_eq_fname);
    append_data_to_file(&spring_constants[i], &eq_bm, 1, bma_eq_fname);
    append_data_to_file(&spring_constants[i], &eq_t , 1, ta_eq_fname);
    append_data_to_file(&spring_constants[i], &eq_um, 1, uma_eq_fname);
  }

  write_config_file(config_eq_fname, CONFIG_OMIT_C, NULL);
  return 0;
}
