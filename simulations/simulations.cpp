#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../dynein_struct.h"

/** Library for simulation code **/

void store_onebound_PEs_callback(void* dyn, State s, void* job_msg, data_union *job_data, int iteration) {
  assert(s == NEARBOUND);
  assert(iteration <= job_data->ob_data.len);
  
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;

  job_data->ob_data.bb[iteration] = dyn_ob->PE_bba;
  job_data->ob_data.bb[iteration] = dyn_ob->PE_bma;
  job_data->ob_data.bb[iteration] = dyn_ob->PE_ta; 
  job_data->ob_data.bb[iteration] = dyn_ob->PE_uma;

  int max_iteration = *((int*) job_msg);

  if (iteration % 100000 == 0) {
    printf("PE calculation progress: %d / %d, %g%%                \r", iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }
  if (iteration == max_iteration) printf("\n");
}

void print_data_to_file(double* data1, double* data2, int iterations, const char* legend, const char* fname) {
  char* buf = (char*)
    malloc((iterations * (12+5+12+1) + 1 + 1) * sizeof(char) + strlen(legend)); // 5 spaced, each double 12 chars, newline
  int buf_offset = 0;
  sprintf(buf, "%s\n", legend);
  buf_offset += strlen(legend) + 1;

  for (int i = 0; i < iterations; i++) {
    sprintf(&buf[buf_offset], "%+.5e     %+.5e\n", data1[i], data2[i]);
    buf_offset += 30;
  }

  FILE* data_file = fopen(fname, "w");
  fputs(buf, data_file);
  fclose(data_file);
}

void get_onebound_PE_correlation_function(generic_data* tau_data, onebound_data* corr_data, int num_corr_datapoints, int min_tau_iter, int max_tau_iter) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double test_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  onebound_data data;
  data.bb = (double*) malloc(num_corr_datapoints * sizeof(double));
  data.bm = (double*) malloc(num_corr_datapoints * sizeof(double));
  data.t  = (double*) malloc(num_corr_datapoints * sizeof(double));
  data.um = (double*) malloc(num_corr_datapoints * sizeof(double));
  data.len = max_tau_iter;

  data_union data_holder;
  data_holder.ob_data = data;

  simulate(max_tau_iter*dt, RAND_INIT_SEED, NEARBOUND, test_position, store_onebound_PEs_callback, (void*) (void*) &max_tau_iter, &data_holder);

  double* PE_bbas = data.bb;
  double* PE_bmas = data.bm;
  double* PE_tas =  data.t;
  double* PE_umas = data.um;

  double PE_bba_ave = get_average(PE_bbas, num_corr_datapoints);
  double PE_bba_var = get_variance(PE_bbas, num_corr_datapoints);
  double PE_bma_ave = get_average(PE_bmas, num_corr_datapoints);
  double PE_bma_var = get_variance(PE_bmas, num_corr_datapoints);
  double PE_ta_ave = get_average(PE_tas, num_corr_datapoints);
  double PE_ta_var = get_variance(PE_tas, num_corr_datapoints);
  double PE_uma_ave = get_average(PE_umas, num_corr_datapoints);
  double PE_uma_var = get_variance(PE_umas, num_corr_datapoints);

  double d_tau_iter = (max_tau_iter - min_tau_iter) / num_corr_datapoints;

  for (int tau_iter = min_tau_iter; tau_iter < max_tau_iter; tau_iter += d_tau_iter) {
    double bba_correlation_sum = 0;
    double bma_correlation_sum = 0;
    double ta_correlation_sum = 0;
    double uma_correlation_sum = 0;
    
    for (int i = 0; i < max_tau_iter - tau_iter; i++) {
      bba_correlation_sum += (PE_bbas[i] - PE_bba_ave) * (PE_bbas[i+tau_iter] - PE_bba_ave);
      bma_correlation_sum += (PE_bmas[i] - PE_bma_ave) * (PE_bmas[i+tau_iter] - PE_bma_ave);
      ta_correlation_sum += (PE_tas[i] - PE_ta_ave) * (PE_tas[i+tau_iter] - PE_ta_ave);
      uma_correlation_sum += (PE_umas[i] - PE_uma_ave) * (PE_umas[i+tau_iter] - PE_uma_ave);
    }

    double bba_correlation = bba_correlation_sum / PE_bba_var / num_corr_datapoints;
    double bma_correlation = bma_correlation_sum / PE_bma_var / num_corr_datapoints;
    double ta_correlation  = ta_correlation_sum / PE_ta_var / num_corr_datapoints;
    double uma_correlation = uma_correlation_sum / PE_uma_var / num_corr_datapoints;

    int iteration = (tau_iter - min_tau_iter) / d_tau_iter;

    assert(iteration <= corr_data->len);

    tau_data->data[iteration] = tau_iter*dt;
    corr_data->bb[iteration] = bba_correlation;
    corr_data->bm[iteration] = bma_correlation;
    corr_data->t[iteration] = ta_correlation;
    corr_data->um[iteration] = uma_correlation;

    if (iteration % 5 == 0) {
      printf("correlation function progress: %g%%                \r", ((double) iteration)  / num_corr_datapoints * 100);
      fflush(NULL);
    }
  }
}

void get_onebound_equipartition_ratio_per_runtime(generic_data* runtime_data, onebound_data* eq_data, int num_eq_datapoints, int min_runtime_iter, int max_runtime_iter) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  onebound_data data;
  data.bb = (double*) malloc(num_eq_datapoints * sizeof(double));
  data.bm = (double*) malloc(num_eq_datapoints * sizeof(double));
  data.t  = (double*) malloc(num_eq_datapoints * sizeof(double));
  data.um = (double*) malloc(num_eq_datapoints * sizeof(double));
  data.len = num_eq_datapoints;

  data_union data_holder;
  data_holder.ob_data = data;
    
  simulate(max_runtime_iter*dt, RAND_INIT_SEED, NEARBOUND, init_position, store_onebound_PEs_callback, (void*) &max_runtime_iter, &data_holder);

  double* bba_PE_data = data_holder.ob_data.bb;
  double* bma_PE_data = data_holder.ob_data.bm;
  double* ta_PE_data =  data_holder.ob_data.t; 
  double* uma_PE_data = data_holder.ob_data.um;

  int d_runtime_iter = (max_runtime_iter - min_runtime_iter) / num_eq_datapoints;

  for (int runtime_iter = min_runtime_iter; runtime_iter < max_runtime_iter; runtime_iter += d_runtime_iter) {
    double bba_eq_ratio = get_average(bba_PE_data, runtime_iter) / (0.5*kb*T);
    double bma_eq_ratio = get_average(bma_PE_data, runtime_iter) / (0.5*kb*T);
    double  ta_eq_ratio = get_average( ta_PE_data, runtime_iter) / (0.5*kb*T);
    double uma_eq_ratio = get_average(uma_PE_data, runtime_iter) / (0.5*kb*T);

    int iteration = (runtime_iter - min_runtime_iter) / d_runtime_iter;
    assert(iteration <= eq_data->len);

    eq_data->bb[iteration] = bba_eq_ratio;
    eq_data->bm[iteration] = bma_eq_ratio;
    eq_data->t[iteration]  = ta_eq_ratio;
    eq_data->um[iteration] = uma_eq_ratio;
    runtime_data->data[iteration] = runtime_iter * dt;

    if (runtime_iter/d_runtime_iter % 1 == 0) {
      printf("equipartition ratio progress: %g%%                \r", ((double) iteration) / num_eq_datapoints * 100);
      fflush(NULL);
    }
  }
}
