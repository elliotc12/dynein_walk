#include <stdio.h>
#include <string.h>
#include <time.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;
  int iterations = 1e6;

  int max_tau_iter = iterations * 0.4;
  int d_tau_iter = 5e4;

  int num_corr_datapoints = max_tau_iter / d_tau_iter;

  const char* bba_corr_title = "Bound binding domain (bba)";
  const char* bma_corr_title = "Bound motor domain (bma)";
  const char* ta_corr_title =  "Tail domain (ta)";
  const char* uma_corr_title = "Unbound motor domain (uma)";

  const char* bba_corr_fname = "bba_pe_correlation.txt";
  const char* bma_corr_fname = "bma_pe_correlation.txt";
  const char* ta_corr_fname = "ta_pe_correlation.txt";
  const char* uma_corr_fname = "uma_pe_correlation.txt";

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

  get_onebound_PE_correlation_function
    (&corr_time_data, &corr_data, d_tau_iter, iterations, max_tau_iter, seeds, seed_len, run_msg);

  print_data_to_file((double*) corr_time_data.data, corr_data.bb, num_corr_datapoints, bba_corr_title, bba_corr_fname);
  print_data_to_file((double*) corr_time_data.data, corr_data.bm, num_corr_datapoints, bma_corr_title, bma_corr_fname);
  print_data_to_file((double*) corr_time_data.data, corr_data.t,  num_corr_datapoints, ta_corr_title, ta_corr_fname);
  print_data_to_file((double*) corr_time_data.data, corr_data.um, num_corr_datapoints, uma_corr_title, uma_corr_fname);

  free (corr_time_data.data);
  free(corr_data.bb); free(corr_data.bm); free(corr_data.t); free(corr_data.um);
  
  return 0;
}
