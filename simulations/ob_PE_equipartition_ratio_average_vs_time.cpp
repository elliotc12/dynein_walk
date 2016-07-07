#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;

  int max_eq_iter = 1e8;
  int min_eq_iter = 0.1*max_eq_iter;
  int num_eq_datapoints = 100;

  const char* bba_eq_title = "Bound binding domain (bba)";
  const char* bma_eq_title = "Bound motor domain (bma)";
  const char* ta_eq_title =  "Tail domain (ta)";
  const char* uma_eq_title = "Unbound motor domain (uma)";

  const char* bba_eq_fname = "bba_pe_equipartition_ratio.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio.txt";

  generic_data eq_time_data;
  onebound_data eq_data, eq_ratio_ave;

  eq_time_data.data = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_time_data.len = num_eq_datapoints;

  eq_data.bb = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.bm = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.t  = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.um = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_data.len = num_eq_datapoints;

  eq_ratio_ave.bb = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.bm = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.t  = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.um = (double*) malloc(num_eq_datapoints * sizeof(double));
  eq_ratio_ave.len = num_eq_datapoints;

  for (int i = 0; i < num_eq_datapoints; i++) {
    eq_ratio_ave.bb[i] = 0;
    eq_ratio_ave.bm[i] = 0;
    eq_ratio_ave.t[i] = 0;
    eq_ratio_ave.um[i] = 0;
  }

  //const int seeds[] = {0, 1, 2, 3, 4};
  const int seeds[] = {0};
  int seed_len = sizeof(seeds) / sizeof(int);

  for (int r = 0; r < seed_len; r++) {
    RAND_INIT_SEED = seeds[r];

    get_onebound_equipartition_ratio_average_per_runtime(&eq_time_data, &eq_data, num_eq_datapoints, min_eq_iter, max_eq_iter);

    for (int i = 0; i < num_eq_datapoints; i++) {
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

  print_data_to_file(eq_time_data.data, eq_ratio_ave.bb, num_eq_datapoints, bba_eq_title, bba_eq_fname);
  print_data_to_file(eq_time_data.data, eq_ratio_ave.bm, num_eq_datapoints, bma_eq_title, bma_eq_fname);
  print_data_to_file(eq_time_data.data, eq_ratio_ave.t, num_eq_datapoints, ta_eq_title, ta_eq_fname);
  print_data_to_file(eq_time_data.data, eq_ratio_ave.um, num_eq_datapoints, uma_eq_title, uma_eq_fname);
  return 0;
}
