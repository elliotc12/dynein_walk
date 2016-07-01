#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../dynein_struct.h"

/** Library for simulation code **/

void store_PEs(void* dyn, State s, void* job_msg, void* job_data, int iteration) {
  assert(s == NEARBOUND);
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
  ((double*) job_data)[4*iteration + 0] = dyn_ob->PE_bba;
  ((double*) job_data)[4*iteration + 1] = dyn_ob->PE_bma;
  ((double*) job_data)[4*iteration + 2] = dyn_ob->PE_ta;
  ((double*) job_data)[4*iteration + 3] = dyn_ob->PE_uma;

  int max_iteration = *((int*) job_msg);

  if (iteration % 100000 == 0) {
    printf("PE calculation progress: %d / %d, %g%%                \r", iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }
  if (iteration == max_iteration) printf("\n");
}


void write_onebound_PE_correlation_function(int iterations, int d_iter, int max_tau_iter) {
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

  int num_correlations = floor(max_tau_iter/d_iter);

  const char* bba_correlation_legend = "PE_nba correlation function";
  const char* bma_correlation_legend = "PE_nma correlation function";
  const char* ta_correlation_legend = "PE_ta correlation function";
  const char* uma_correlation_legend = "PE_uma correlation function";
  
  char* bba_correlation_buf = (char*)
    malloc((num_correlations * (12+5+12+1) + 1 + 1) * sizeof(char) + strlen(bba_correlation_legend));
  char* bma_correlation_buf = (char*)
    malloc((num_correlations * (12+5+12+1) + 1 + 1) * sizeof(char) + strlen(bma_correlation_legend));
  char* ta_correlation_buf = (char*)
    malloc((num_correlations * (12+5+12+1) + 1 + 1) * sizeof(char) + strlen(ta_correlation_legend));
  char* uma_correlation_buf = (char*)
    malloc((num_correlations * (12+5+12+1) + 1 + 1) * sizeof(char) + strlen(uma_correlation_legend));

  int bba_correlation_buf_offset = 0;
  int bma_correlation_buf_offset = 0;
  int ta_correlation_buf_offset = 0;
  int uma_correlation_buf_offset = 0;
  
  sprintf(bba_correlation_buf, "%s\n", bba_correlation_legend);
  sprintf(bma_correlation_buf, "%s\n", bma_correlation_legend);
  sprintf(ta_correlation_buf, "%s\n", ta_correlation_legend);
  sprintf(uma_correlation_buf, "%s\n", uma_correlation_legend);
  
  bba_correlation_buf_offset += strlen(bba_correlation_legend) + 1;
  bma_correlation_buf_offset += strlen(bma_correlation_legend) + 1;
  ta_correlation_buf_offset += strlen(ta_correlation_legend) + 1;
  uma_correlation_buf_offset += strlen(uma_correlation_legend) + 1;

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

    sprintf(&bba_correlation_buf[bba_correlation_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, bba_correlation);
    sprintf(&bma_correlation_buf[bma_correlation_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, bma_correlation);
    sprintf(&ta_correlation_buf[ta_correlation_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, ta_correlation);
    sprintf(&uma_correlation_buf[uma_correlation_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, uma_correlation);

    bba_correlation_buf_offset += 30;
    bma_correlation_buf_offset += 30;
    ta_correlation_buf_offset += 30;
    uma_correlation_buf_offset += 30;

    if (tau_iter/d_iter % 5 == 0) {
    printf("correlation function progress: %d / %d, %g%%                \r",
	   tau_iter, max_tau_iter, ((double) tau_iter) / max_tau_iter * 100);
    fflush(NULL);
    }
  }
  
  FILE* bba_correlation_data_file = fopen("pe_bba_correlation_function.txt", "w");
  FILE* bma_correlation_data_file = fopen("pe_bma_correlation_function.txt", "w");
  FILE* ta_correlation_data_file = fopen("pe_ta_correlation_function.txt", "w");
  FILE* uma_correlation_data_file = fopen("pe_uma_correlation_function.txt", "w");
  
  fputs(bba_correlation_buf, bba_correlation_data_file);
  fputs(bma_correlation_buf, bma_correlation_data_file);
  fputs(ta_correlation_buf, ta_correlation_data_file);
  fputs(uma_correlation_buf, uma_correlation_data_file);
  
  fclose(bba_correlation_data_file);
  fclose(bma_correlation_data_file);
  fclose(ta_correlation_data_file);
  fclose(uma_correlation_data_file);
}

void write_onebound_equipartition_ratio_per_tau(int iterations, int d_iter, int min_tau_iter, int max_tau_iter) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  
  double runtime = dt*iterations;
  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double test_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};
  void* data = malloc(iterations * sizeof(double) * 4);
    
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, test_position, store_PEs, (void*) &iterations, data);

  double* bba_PEs = (double*) malloc(iterations * sizeof(double));
  double* bma_PEs = (double*) malloc(iterations * sizeof(double));
  double*  ta_PEs = (double*) malloc(iterations * sizeof(double));
  double* uma_PEs = (double*) malloc(iterations * sizeof(double));

  for (int i = 0; i < iterations; i++) {
    bba_PEs[i] = ((double*) data)[4*i + 0];
    bma_PEs[i] = ((double*) data)[4*i + 1];
    ta_PEs[i]  = ((double*) data)[4*i + 2];
    uma_PEs[i] = ((double*) data)[4*i + 3];
  }

  const char* bba_eq_ratio_legend = "PE_nba / 0.5*kb*T ratio";
  const char* bma_eq_ratio_legend = "PE_nma / 0.5*kb*T ratio";
  const char* ta_eq_ratio_legend = "PE_ta / 0.5*kb*T ratio";
  const char* uma_eq_ratio_legend = "PE_uma / 0.5*kb*T ratio";

  int num_iterations = floor(max_tau_iter/d_iter);

  char* bba_eq_ratio_buf = (char*)
    malloc((num_iterations * 30 + 1 + 1) * sizeof(char) + strlen(bba_eq_ratio_legend));
  char* bma_eq_ratio_buf = (char*)
    malloc((num_iterations * 30 + 1 + 1) * sizeof(char) + strlen(bma_eq_ratio_legend));
  char* ta_eq_ratio_buf = (char*)
    malloc((num_iterations * 30 + 1 + 1) * sizeof(char) + strlen(ta_eq_ratio_legend));
  char* uma_eq_ratio_buf = (char*)
    malloc((num_iterations * 30 + 1 + 1) * sizeof(char) + strlen(uma_eq_ratio_legend));

  int bba_eq_ratio_buf_offset = 0;
  int bma_eq_ratio_buf_offset = 0;
  int ta_eq_ratio_buf_offset = 0;
  int uma_eq_ratio_buf_offset = 0;

  sprintf(bba_eq_ratio_buf, "%s\n", bba_eq_ratio_legend);
  sprintf(bma_eq_ratio_buf, "%s\n", bma_eq_ratio_legend);
  sprintf(ta_eq_ratio_buf, "%s\n", ta_eq_ratio_legend);
  sprintf(uma_eq_ratio_buf, "%s\n", uma_eq_ratio_legend);
  
  bba_eq_ratio_buf_offset += strlen(bba_eq_ratio_legend) + 1;
  bma_eq_ratio_buf_offset += strlen(bma_eq_ratio_legend) + 1;
  ta_eq_ratio_buf_offset += strlen(ta_eq_ratio_legend) + 1;
  uma_eq_ratio_buf_offset += strlen(uma_eq_ratio_legend) + 1;

  for (int tau_iter=min_tau_iter; tau_iter < max_tau_iter; tau_iter += d_iter) {
    double bba_eq_ratio = get_average(bba_PEs, tau_iter) / (0.5*kb*T);
    double bma_eq_ratio = get_average(bma_PEs, tau_iter) / (0.5*kb*T);
    double  ta_eq_ratio = get_average( ta_PEs, tau_iter) / (0.5*kb*T);
    double uma_eq_ratio = get_average(uma_PEs, tau_iter) / (0.5*kb*T);
    
    sprintf(&bba_eq_ratio_buf[bba_eq_ratio_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, bba_eq_ratio);
    sprintf(&bma_eq_ratio_buf[bma_eq_ratio_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, bma_eq_ratio);
    sprintf(&ta_eq_ratio_buf[ta_eq_ratio_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, ta_eq_ratio);
    sprintf(&uma_eq_ratio_buf[uma_eq_ratio_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, uma_eq_ratio);
    
    bba_eq_ratio_buf_offset += 30;
    bma_eq_ratio_buf_offset += 30;
    ta_eq_ratio_buf_offset += 30;
    uma_eq_ratio_buf_offset += 30;

    if (tau_iter/d_iter % 5 == 0) {
      printf("equipartition ratio progress: %d / %d, %g%%                \r",
	     tau_iter, max_tau_iter, ((double) tau_iter) / max_tau_iter * 100);
      fflush(NULL);
    }
  }

  FILE* bba_eq_ratio_data_file = fopen("bba_equipartition_ratio.txt", "w");
  FILE* bma_eq_ratio_data_file = fopen("bma_equipartition_ratio.txt", "w");
  FILE* ta_eq_ratio_data_file = fopen("ta_equipartition_ratio.txt", "w");
  FILE* uma_eq_ratio_data_file = fopen("uma_equipartition_ratio.txt", "w");
  
  fputs(bba_eq_ratio_buf, bba_eq_ratio_data_file);
  fputs(bma_eq_ratio_buf, bma_eq_ratio_data_file);
  fputs(ta_eq_ratio_buf, ta_eq_ratio_data_file);
  fputs(uma_eq_ratio_buf, uma_eq_ratio_data_file);
  
  fclose(bba_eq_ratio_data_file);
  fclose(bma_eq_ratio_data_file);
  fclose(ta_eq_ratio_data_file);
  fclose(uma_eq_ratio_data_file);
}
