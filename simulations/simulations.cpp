#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../dynein_struct.h"

/** Library for simulation code **/


// David -- naming convention for callback fn-ish things like this?
void store_PEs(void* dyn, State s, void* job_msg, void* job_data, int iteration) {
  assert(s == NEARBOUND);
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
  ((double*) job_data)[4*iteration + 0] = dyn_ob->PE_bba;
  ((double*) job_data)[4*iteration + 1] = dyn_ob->PE_bma;
  ((double*) job_data)[4*iteration + 2] = dyn_ob->PE_ta;
  ((double*) job_data)[4*iteration + 3] = dyn_ob->PE_uma;

  int max_iteration = *((int*) job_msg);

  if (iteration % 100000 == 0) {
    //printf("PE calculation progress: %d / %d, %g%%                \r", iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    printf("PE calculation progress: %d / %d, %g%% rand: %p         \n", iteration, max_iteration, ((double) iteration) / max_iteration * 100, dyn_ob->rand);
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

void get_onebound_PE_correlation_function(double* time_data, double* corr_data, int iterations, int d_iter, int max_tau_iter) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  
  double runtime = dt*iterations;
  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double test_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};
  void* data = malloc(4 * iterations * sizeof(double));
  
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, test_position, store_PEs, (void*) &iterations, data);

  double* PE_bbas = (double*) malloc(iterations * sizeof(double));
  double* PE_bmas = (double*) malloc(iterations * sizeof(double));
  double* PE_tas = (double*) malloc(iterations * sizeof(double));
  double* PE_umas = (double*) malloc(iterations * sizeof(double));
  
  for (int i = 0; i < iterations; i++) {
    PE_bbas[i] = ((double*) data)[4*i + 0];
    PE_bmas[i] = ((double*) data)[4*i + 1];
    PE_tas[i]  = ((double*) data)[4*i + 2];
    PE_umas[i] = ((double*) data)[4*i + 3];
  }

  double PE_bba_ave = get_average(PE_bbas, iterations);
  double PE_bba_var = get_variance(PE_bbas, iterations);
  double PE_bma_ave = get_average(PE_bmas, iterations);
  double PE_bma_var = get_variance(PE_bmas, iterations);
  double PE_ta_ave = get_average(PE_tas, iterations);
  double PE_ta_var = get_variance(PE_tas, iterations);
  double PE_uma_ave = get_average(PE_umas, iterations);
  double PE_uma_var = get_variance(PE_umas, iterations);

  for (int tau_iter=0; tau_iter < max_tau_iter; tau_iter += d_iter) {
    double bba_correlation_sum = 0;
    double bma_correlation_sum = 0;
    double ta_correlation_sum = 0;
    double uma_correlation_sum = 0;
    
    for (int i = 0; i < iterations - tau_iter; i++) {
      bba_correlation_sum += (PE_bbas[i] - PE_bba_ave) * (PE_bbas[i+tau_iter] - PE_bba_ave);
      bma_correlation_sum += (PE_bmas[i] - PE_bma_ave) * (PE_bmas[i+tau_iter] - PE_bma_ave);
      ta_correlation_sum += (PE_tas[i] - PE_ta_ave) * (PE_tas[i+tau_iter] - PE_ta_ave);
      uma_correlation_sum += (PE_umas[i] - PE_uma_ave) * (PE_umas[i+tau_iter] - PE_uma_ave);
    }
    
    int num_iterations = iterations - tau_iter;

    double bba_correlation = bba_correlation_sum / PE_bba_var / num_iterations;
    double bma_correlation = bma_correlation_sum / PE_bma_var / num_iterations;
    double ta_correlation  = ta_correlation_sum / PE_ta_var / num_iterations;
    double uma_correlation = uma_correlation_sum / PE_uma_var / num_iterations;

    time_data[tau_iter / d_iter] = tau_iter*dt;
    corr_data[4*(tau_iter / d_iter) + 0] = bba_correlation;
    corr_data[4*(tau_iter / d_iter) + 1] = bma_correlation;
    corr_data[4*(tau_iter / d_iter) + 2] = ta_correlation;
    corr_data[4*(tau_iter / d_iter) + 3] = uma_correlation;

    if (tau_iter/d_iter % 5 == 0) {
    printf("correlation function progress: %d / %d, %g%%                \r",
	   tau_iter, max_tau_iter, ((double) tau_iter) / max_tau_iter * 100);
    fflush(NULL);
    }
  }
}

void get_onebound_equipartition_ratio_per_runtime(double* runtime_data, double* eq_data, int d_runtime_iter, int min_runtime_iter, int max_runtime_iter) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  
  double runtime = dt*max_runtime_iter;
  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};
  void* data = malloc(ceil(max_runtime_iter/d_runtime_iter) * sizeof(double) * 4);
    
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, init_position, store_PEs, (void*) &max_runtime_iter, data);

  int num_PE_datapoints = floor((max_runtime_iter - min_runtime_iter) / d_runtime_iter);

  double* bba_PE_data = (double*) malloc(num_PE_datapoints * sizeof(double));
  double* bma_PE_data = (double*) malloc(num_PE_datapoints * sizeof(double));
  double* ta_PE_data =  (double*) malloc(num_PE_datapoints * sizeof(double));
  double* uma_PE_data = (double*) malloc(num_PE_datapoints * sizeof(double));

  for (int i = 0; i < num_PE_datapoints; i++) {          // unpack equipartition data
    bba_PE_data[i] = ((double*) data)[4*i + 0];
    bma_PE_data[i] = ((double*) data)[4*i + 1];
    ta_PE_data[i]  = ((double*) data)[4*i + 2];
    uma_PE_data[i] = ((double*) data)[4*i + 3];
  }

  for (int runtime_iter = min_runtime_iter; runtime_iter < max_runtime_iter; runtime_iter += d_runtime_iter) {
    double bba_eq_ratio = get_average(bba_PE_data, runtime_iter) / (0.5*kb*T);
    double bma_eq_ratio = get_average(bma_PE_data, runtime_iter) / (0.5*kb*T);
    double  ta_eq_ratio = get_average( ta_PE_data, runtime_iter) / (0.5*kb*T);
    double uma_eq_ratio = get_average(uma_PE_data, runtime_iter) / (0.5*kb*T);

    eq_data[4*((runtime_iter - min_runtime_iter) / d_runtime_iter) + 0] = bba_eq_ratio;
    eq_data[4*((runtime_iter - min_runtime_iter) / d_runtime_iter) + 1] = bma_eq_ratio;
    eq_data[4*((runtime_iter - min_runtime_iter) / d_runtime_iter) + 2] = ta_eq_ratio;
    eq_data[4*((runtime_iter - min_runtime_iter) / d_runtime_iter) + 3] = uma_eq_ratio;

    runtime_data[(runtime_iter - min_runtime_iter) / d_runtime_iter] = runtime_iter * dt;

    if (runtime_iter/d_runtime_iter % 1 == 0) {
      //printf("equipartition ratio progress: %d / %d, %g%%                \r",
      printf("equipartition ratio progress: %d / %d, %g%%                \n",
	     runtime_iter - min_runtime_iter, max_runtime_iter - min_runtime_iter, ((double) (runtime_iter - min_runtime_iter)) / (max_runtime_iter - min_runtime_iter) * 100);
      fflush(NULL);
    }
  }
}
