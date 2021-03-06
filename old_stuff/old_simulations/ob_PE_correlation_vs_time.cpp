#include <cassert>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main(int argc, char** argv) {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 1000;
  int iterations = 1e6;

  int max_tau_iter = iterations * 0.4;
  int num_corr_datapoints = 1000;

  const char* bba_corr_title = "Bound binding";
  const char* bma_corr_title = "Bound motor";
  const char* ta_corr_title =  "Tail domain";
  const char* uma_corr_title = "Unbound motor";

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }
  char* f_appended_name = argv[1];
  char bba_corr_fname[200];
  char bma_corr_fname[200];
  char ta_corr_fname[200];
  char uma_corr_fname[200];
  char config_corr_fname[200];

  strcpy(bba_corr_fname, "data/bba_pe_correlation_");
  strcpy(bma_corr_fname, "data/bma_pe_correlation_");
  strcpy(ta_corr_fname,  "data/ta_pe_correlation_");
  strcpy(uma_corr_fname, "data/uma_pe_correlation_");
  strcpy(uma_corr_fname, "data/config_pe_correlation_");

  strcat(bba_corr_fname, f_appended_name);
  strcat(bma_corr_fname, f_appended_name);
  strcat(ta_corr_fname, f_appended_name);
  strcat(uma_corr_fname, f_appended_name);
  strcat(config_corr_fname, f_appended_name);

  strcat(bba_corr_fname, ".txt");
  strcat(bma_corr_fname, ".txt");
  strcat(ta_corr_fname, ".txt");
  strcat(uma_corr_fname, ".txt");
  strcat(config_corr_fname, ".txt");

  prepare_data_file(bba_corr_title, bba_corr_fname);
  prepare_data_file(bma_corr_title, bma_corr_fname);
  prepare_data_file(ta_corr_title,  ta_corr_fname);
  prepare_data_file(uma_corr_title, uma_corr_fname);

  generic_data corr_time_data;
  onebound_data corr_data;

  corr_time_data.data = (void*) malloc(num_corr_datapoints * sizeof(double));
  corr_time_data.len = num_corr_datapoints;

  corr_data.bb = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.bm = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.t  = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.um = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.len = num_corr_datapoints;

  const int seeds[] = {0};
  int seed_len = sizeof(seeds) / sizeof(int);

  char run_msg[512];
  const char* run_msg_base = "correlation calc (";
  strcpy(run_msg, run_msg_base);

  get_onebound_PE_correlation_function(&corr_time_data, &corr_data, num_corr_datapoints,
     iterations,max_tau_iter,seeds, seed_len, run_msg);

  append_data_to_file((double*)corr_time_data.data,corr_data.bb, num_corr_datapoints, bba_corr_fname);
  append_data_to_file((double*)corr_time_data.data,corr_data.bm, num_corr_datapoints, bma_corr_fname);
  append_data_to_file((double*)corr_time_data.data,corr_data.t , num_corr_datapoints, ta_corr_fname);
  append_data_to_file((double*)corr_time_data.data,corr_data.um, num_corr_datapoints, uma_corr_fname);
  write_config_file(config_corr_fname, NULL, NULL);

  free (corr_time_data.data);
  free(corr_data.bb); free(corr_data.bm); free(corr_data.t); free(corr_data.um);
  
  return 0;
}
