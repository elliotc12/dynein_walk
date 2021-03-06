#include <cassert>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main(int argc, char** argv) {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 1000;
  long long max_eq_iter = 1e10;
  long long min_eq_iter = 1e2;
  int num_eq_datapoints = 500;

  Lt = 15;
  Ls = 15;
  fake_radius_t = 1.5;
  fake_radius_m = 1.5;
  fake_radius_b = 1.5;
  gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
  gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
  gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

  cb = 0.5;
  cm = 0.5;
  ct = 0.5;

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

  strcpy(bba_eq_fname, "data/bba_pe_log_equipartition_vs_log_iterations_");
  strcpy(bma_eq_fname, "data/bma_pe_log_equipartition_vs_log_iterations_");
  strcpy(ta_eq_fname,  "data/ta_pe_log_equipartition_vs_log_iterations_");
  strcpy(config_eq_fname, "data/config_pe_log_equipartition_vs_log_iterations_");
  

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

  generic_data eq_time_data;
  onebound_data eq_data;

  eq_time_data.data = (void*) malloc(num_eq_datapoints * sizeof(double));
  eq_time_data.len = num_eq_datapoints;

  eq_data.bb = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.bm = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.t  = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.um = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.len = num_eq_datapoints;

  const int seeds[] = {0};
  int seed_len = sizeof(seeds) / sizeof(int);

  char run_msg[512];
  const char* run_msg_base = "ET v iter (";
  strcpy(run_msg, run_msg_base);

  get_onebound_equipartition_ratio_average_per_runtime(&eq_time_data, &eq_data, num_eq_datapoints,
        min_eq_iter, max_eq_iter, seeds, seed_len, run_msg);

  for (int i = 0; i < num_eq_datapoints; i++) {
    eq_data.bb[i] = std::abs(eq_data.bb[i] - 1);
    eq_data.bm[i] = std::abs(eq_data.bm[i] - 1);
    eq_data.t[i]  = std::abs(eq_data.t[i] - 1);
    eq_data.um[i] = std::abs(eq_data.um[i] - 1);
    ((double*) eq_time_data.data)[i] = ((double*) eq_time_data.data)[i]/dt;
  }

  append_data_to_file((double*) eq_time_data.data, eq_data.bb, num_eq_datapoints, bba_eq_fname);
  append_data_to_file((double*) eq_time_data.data, eq_data.bm, num_eq_datapoints, bma_eq_fname);
  append_data_to_file((double*) eq_time_data.data, eq_data.t , num_eq_datapoints, ta_eq_fname);
  append_data_to_file((double*) eq_time_data.data, eq_data.um, num_eq_datapoints, uma_eq_fname);
  write_config_file(config_eq_fname, NULL, NULL);

  free(eq_time_data.data);
  free(eq_data.bb); free(eq_data.bm); free(eq_data.t); free(eq_data.um);
  
  return 0;
}
