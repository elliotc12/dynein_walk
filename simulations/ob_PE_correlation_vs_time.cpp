#include <stdio.h>
#include <time.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;
  int iterations = 1e7;

  int max_tau_iter = iterations * 0.4;
  int num_corr_datapoints = 50;

  const char* bba_corr_title = "Bound binding domain (bba)";
  const char* bma_corr_title = "Bound motor domain (bma)";
  const char* ta_corr_title =  "Tail domain (ta)";
  const char* uma_corr_title = "Unbound motor domain (uma)";

  const char* bba_corr_fname = "bba_pe_correlation.txt";
  const char* bma_corr_fname = "bma_pe_correlation.txt";
  const char* ta_corr_fname = "ta_pe_correlation.txt";
  const char* uma_corr_fname = "uma_pe_correlation.txt";

  generic_data corr_time_data;
  onebound_data corr_data, corr_ave;

  corr_time_data.data = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_time_data.len = num_corr_datapoints;
  
  corr_data.bb = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.bm = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.t  = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.um = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_data.len = num_corr_datapoints;

  corr_ave.bb = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.bm = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.t  = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.um = (double*) malloc(num_corr_datapoints * sizeof(double));
  corr_ave.len = num_corr_datapoints;

  for (int i = 0; i < num_corr_datapoints; i++) {
    corr_ave.bb[i] = 0;
    corr_ave.bm[i] = 0;
    corr_ave.t[i] = 0;
    corr_ave.um[i] = 0;
  }

  //const int seeds[] = {0, 1, 2, 3, 4};
  const int seeds[] = {0};
  int seed_len = sizeof(seeds) / sizeof(int);

  for (int r = 0; r < seed_len; r++) {
    RAND_INIT_SEED = seeds[r];

    get_onebound_PE_correlation_function(&corr_time_data, &corr_data, num_corr_datapoints, iterations, max_tau_iter);

    for (int i = 0; i < num_corr_datapoints; i++) {
      corr_ave.bb[i] += corr_data.bb[i];
      corr_ave.bm[i] += corr_data.bm[i];
      corr_ave.t[i]  += corr_data.t[i];
      corr_ave.um[i] += corr_data.um[i];
    }
  }
  
  for (int i = 0; i < num_corr_datapoints; i++) {
    corr_ave.bb[i] /= seed_len;
    corr_ave.bm[i] /= seed_len;
    corr_ave.t[i]  /= seed_len;
    corr_ave.um[i] /= seed_len;
  }

  print_data_to_file(corr_time_data.data, corr_ave.bb, num_corr_datapoints, bba_corr_title, bba_corr_fname);
  print_data_to_file(corr_time_data.data, corr_ave.bm, num_corr_datapoints, bma_corr_title, bma_corr_fname);
  print_data_to_file(corr_time_data.data, corr_ave.t,  num_corr_datapoints, ta_corr_title, ta_corr_fname);
  print_data_to_file(corr_time_data.data, corr_ave.um, num_corr_datapoints, uma_corr_title, uma_corr_fname);
  return 0;
}
