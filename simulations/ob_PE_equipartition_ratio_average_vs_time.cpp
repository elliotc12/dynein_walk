#include <stdio.h>
#include <string.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 1000;
  int max_eq_iter = 1e6;
  int min_eq_iter = 0.1*max_eq_iter;
  int d_runtime_iter = 1e5;
  
  int num_eq_datapoints = (max_eq_iter - min_eq_iter) / d_runtime_iter;

  const char* bba_eq_title = "Bound binding";
  const char* bma_eq_title = "Bound motor";
  const char* ta_eq_title =  "Tail domain";
  const char* uma_eq_title = "Unbound motor";

  const char* bba_eq_fname = "bba_pe_equipartition_ratio_v_time.txt";
  const char* bma_eq_fname = "bma_pe_equipartition_ratio_v_time.txt";
  const char* ta_eq_fname =  "ta_pe_equipartition_ratio_v_time.txt";
  const char* uma_eq_fname = "uma_pe_equipartition_ratio_v_time.txt";

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
  const char* run_msg_base = "timecalc (";
  strcpy(run_msg, run_msg_base);

  get_onebound_equipartition_ratio_average_per_runtime(&eq_time_data, &eq_data, d_runtime_iter, min_eq_iter, max_eq_iter, seeds, seed_len, run_msg);

  append_data_to_file((double*) eq_time_data.data,eq_data.bb,num_eq_datapoints, bba_eq_fname);
  append_data_to_file((double*) eq_time_data.data,eq_data.bm,num_eq_datapoints, bma_eq_fname);
  append_data_to_file((double*) eq_time_data.data,eq_data.t ,num_eq_datapoints, ta_eq_fname);
  append_data_to_file((double*) eq_time_data.data,eq_data.um,num_eq_datapoints, uma_eq_fname);

  free (eq_time_data.data);
  free(eq_data.bb); free(eq_data.bm); free(eq_data.t); free(eq_data.um);
  
  return 0;
}
