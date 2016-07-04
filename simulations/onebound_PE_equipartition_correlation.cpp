#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"
#include "simulations.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;
  int iterations = 1e6;
  int d_iter = 1e4;
  int min_tau_iter = floor(iterations * 0.1);
  int max_tau_iter = floor(iterations * 0.5);

  int min_runtime_iter = 10*min_tau_iter;
  int max_runtime_iter = 10*max_tau_iter;
  int d_runtime_iter = 10*d_iter;
  int num_eq_datapoints = ceil((max_runtime_iter - min_runtime_iter) / d_runtime_iter);
  int num_correlations = ceil((max_tau_iter - min_tau_iter) / d_iter);

  const char* bba_corr_title = "Bound binding domain (bba)";
  const char* bma_corr_title = "Bound motor domain (bma)";
  const char* ta_corr_title =  "Tail domain (ta)";
  const char* uma_corr_title = "Unbound motor domain (uma)";

  const char* bba_corr_fname = "bba_pe_correlation.txt";
  const char* bma_corr_fname = "bma_pe_correlation.txt";
  const char* ta_corr_fname = "ta_pe_correlation.txt";
  const char* uma_corr_fname = "uma_pe_correlation.txt";

  const char* bba_eq_title = "Bound binding domain (bba)";
  const char* bma_eq_title = "Bound motor domain (bma)";
  const char* ta_eq_title =  "Tail domain (ta)";
  const char* uma_eq_title = "Unbound motor domain (uma)";

  const char* bba_eq_fname = "bba_pe_equipartition_ratio.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio.txt";

  double* tau_data = (double*) malloc(num_correlations * sizeof(double));
  double* corr_data = (double*) malloc(4 * num_correlations * sizeof(double));
  
  double* runtime_data = (double*) malloc(num_eq_datapoints * sizeof(double));
  double* eq_data = (double*) malloc(4 * num_eq_datapoints * sizeof(double));

  double* bba_eq_ratio_ave = (double*) malloc(num_eq_datapoints * sizeof(double));
  double* bma_eq_ratio_ave = (double*) malloc(num_eq_datapoints * sizeof(double));
  double* ta_eq_ratio_ave =  (double*) malloc(num_eq_datapoints * sizeof(double));
  double* uma_eq_ratio_ave = (double*) malloc(num_eq_datapoints * sizeof(double));

  double* bba_corr_ave = (double*) malloc(num_correlations * sizeof(double));
  double* bma_corr_ave = (double*) malloc(num_correlations * sizeof(double));
  double* ta_corr_ave =  (double*) malloc(num_correlations * sizeof(double));
  double* uma_corr_ave = (double*) malloc(num_correlations * sizeof(double));

  for (int i = 0; i < num_eq_datapoints; i++) {
    bba_eq_ratio_ave[i] = 0;
    bma_eq_ratio_ave[i] = 0;
    ta_eq_ratio_ave[i] = 0;
    uma_eq_ratio_ave[i] = 0;
  }

  for (int i = 0; i < num_correlations; i++) {
    bba_corr_ave[i] = 0;
    bma_corr_ave[i] = 0;
    ta_corr_ave[i] = 0;
    uma_corr_ave[i] = 0;
  }

  const int seeds[] = {0, 1, 2, 3, 4};
  int seed_len = sizeof(seeds) / sizeof(int);

  for (int r = 0; r < seed_len; r++) {
    RAND_INIT_SEED = seeds[r];

    get_onebound_PE_correlation_function(tau_data, corr_data, iterations, d_iter, max_tau_iter);
    get_onebound_equipartition_ratio_per_runtime(runtime_data, eq_data, d_runtime_iter, min_runtime_iter, max_runtime_iter);

    for (int i = 0; i < num_correlations; i++) {          // unpack correlation data
      bba_corr_ave[i] += ((double*) corr_data)[4*i + 0];
      bma_corr_ave[i] += ((double*) corr_data)[4*i + 1];
      ta_corr_ave[i]  += ((double*) corr_data)[4*i + 2];
      uma_corr_ave[i] += ((double*) corr_data)[4*i + 3];
    }

    for (int i = 0; i < num_eq_datapoints; i++) {          // unpack equipartition data
      bba_eq_ratio_ave[i] += ((double*) eq_data)[4*i + 0];
      bma_eq_ratio_ave[i] += ((double*) eq_data)[4*i + 1];
      ta_eq_ratio_ave[i]  += ((double*) eq_data)[4*i + 2];
      uma_eq_ratio_ave[i] += ((double*) eq_data)[4*i + 3];
    }
  }

  for (int i = 0; i < num_eq_datapoints; i++) {
    bba_eq_ratio_ave[i] /= seed_len;
    bma_eq_ratio_ave[i] /= seed_len;
    ta_eq_ratio_ave[i] /= seed_len;
    uma_eq_ratio_ave[i] /= seed_len;
  }
  
  for (int i = 0; i < num_correlations; i++) {
    bba_corr_ave[i] /= seed_len;
    bma_corr_ave[i] /= seed_len;
    ta_corr_ave[i] /= seed_len;
    uma_corr_ave[i] /= seed_len;
  }
  

  print_data_to_file(tau_data, bba_corr_ave, num_correlations, bba_corr_title, bba_corr_fname);
  print_data_to_file(tau_data, bma_corr_ave, num_correlations, bma_corr_title, bma_corr_fname);
  print_data_to_file(tau_data, ta_corr_ave, num_correlations, ta_corr_title, ta_corr_fname);
  print_data_to_file(tau_data, uma_corr_ave, num_correlations, uma_corr_title, uma_corr_fname);

  print_data_to_file(runtime_data, bba_eq_ratio_ave, num_eq_datapoints, bba_eq_title, bba_eq_fname);
  print_data_to_file(runtime_data, bma_eq_ratio_ave, num_eq_datapoints, bma_eq_title, bma_eq_fname);
  print_data_to_file(runtime_data, ta_eq_ratio_ave, num_eq_datapoints, ta_eq_title, ta_eq_fname);
  print_data_to_file(runtime_data, uma_eq_ratio_ave, num_eq_datapoints, uma_eq_title, uma_eq_fname);
  return 0;
}
