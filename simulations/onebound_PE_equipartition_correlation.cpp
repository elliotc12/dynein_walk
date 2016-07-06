#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;
  int iterations = 1e6;

  int min_tau_iter = iterations * 0.1;
  int max_tau_iter = iterations * 0.5;
  int num_corr_datapoints = 50;

  int min_runtime_iter = 10*min_tau_iter;
  int max_runtime_iter = 10*max_tau_iter;
  int num_eq_datapoints = 100;

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

  generic_data tau_data, runtime_data;
  onebound_data corr_data, eq_data, eq_ratio_ave, corr_ave;

  tau_data.data = (double*) malloc(num_corr_datapoints * sizeof(double));
  runtime_data.data = (double*) malloc(num_eq_datapoints * sizeof(double));
  
  corr_data.bb = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.bm = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.t  = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.um = (double*) malloc(num_corr_datapoints * sizeof(double));
  
  eq_data.bb = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.bm = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.t  = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.um = (double*) malloc(num_eq_datapoints * sizeof(double));

  eq_ratio_ave.bb = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.bm = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.t  = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.um = (double*) malloc(num_eq_datapoints * sizeof(double));

  corr_ave.bb = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.bm = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.t  = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.um = (double*) malloc(num_corr_datapoints * sizeof(double));

  for (int i = 0; i < num_eq_datapoints; i++) {
    eq_ratio_ave.bb[i] = 0;
    eq_ratio_ave.bm[i] = 0;
    eq_ratio_ave.t[i] = 0;
    eq_ratio_ave.um[i] = 0;
  }

  for (int i = 0; i < num_corr_datapoints; i++) {
    corr_ave.bb[i] = 0;
    corr_ave.bm[i] = 0;
    corr_ave.t[i] = 0;
    corr_ave.um[i] = 0;
  }

  const int seeds[] = {0, 1, 2, 3, 4};
  int seed_len = sizeof(seeds) / sizeof(int);

  for (int r = 0; r < seed_len; r++) {
    RAND_INIT_SEED = seeds[r];

    get_onebound_PE_correlation_function(&tau_data, &corr_data, num_corr_datapoints, min_tau_iter, max_tau_iter);
    get_onebound_equipartition_ratio_per_runtime(&runtime_data, &eq_data, num_eq_datapoints, min_runtime_iter, max_runtime_iter);

    for (int i = 0; i < num_corr_datapoints; i++) {          // unpack correlation data
      corr_ave.bb[i] += corr_data.bb[i];
      corr_ave.bm[i] += corr_data.bm[i];
      corr_ave.t[i]  += corr_data.t[i];
      corr_ave.um[i] += corr_data.um[i];
    }

    for (int i = 0; i < num_eq_datapoints; i++) {         // unpack equipartition data
      eq_ratio_ave.bb[i] += eq_data.bb[i];
      eq_ratio_ave.bm[i] += eq_data.bm[i];
      eq_ratio_ave.t[i]  += eq_data.t[i];
      eq_ratio_ave.um[i] += eq_data.um[i];
    }
  }

  for (int i = 0; i < num_eq_datapoints; i++) {
    eq_ratio_ave.bb[i] /= seed_len;
    eq_ratio_ave.bm[i] /= seed_len;
    eq_ratio_ave.t[i]  /= seed_len;
    eq_ratio_ave.um[i] /= seed_len;
  }
  
  for (int i = 0; i < num_corr_datapoints; i++) {
    corr_ave.bb[i] /= seed_len;
    corr_ave.bm[i] /= seed_len;
    corr_ave.t[i]  /= seed_len;
    corr_ave.um[i] /= seed_len;
  }

  print_data_to_file(tau_data.data, corr_ave.bb, num_corr_datapoints, bba_corr_title, bba_corr_fname);
  print_data_to_file(tau_data.data, corr_ave.bm, num_corr_datapoints, bma_corr_title, bma_corr_fname);
  print_data_to_file(tau_data.data, corr_ave.t,  num_corr_datapoints, ta_corr_title, ta_corr_fname);
  print_data_to_file(tau_data.data, corr_ave.um, num_corr_datapoints, uma_corr_title, uma_corr_fname);

  print_data_to_file(runtime_data.data, eq_ratio_ave.bb, num_eq_datapoints, bba_eq_title, bba_eq_fname);
  print_data_to_file(runtime_data.data, eq_ratio_ave.bm, num_eq_datapoints, bma_eq_title, bma_eq_fname);
  print_data_to_file(runtime_data.data, eq_ratio_ave.t, num_eq_datapoints, ta_eq_title, ta_eq_fname);
  print_data_to_file(runtime_data.data, eq_ratio_ave.um, num_eq_datapoints, uma_eq_title, uma_eq_fname);
  return 0;
}
