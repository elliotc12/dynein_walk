#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../dynein_struct.h"

/** Library for simulation code **/

void store_onebound_PEs_callback(void* dyn, State s, void* job_msg, data_union *job_data, long long iteration) {
  assert(s == NEARBOUND);
  assert(iteration < job_data->ob_data.len);
  long long max_iteration = (long long) ((double*) job_msg)[0];
  double start_time= ((double*) job_msg)[1];
  
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
  job_data->ob_data.bb[iteration] = dyn_ob->PE_bba;
  job_data->ob_data.bm[iteration] = dyn_ob->PE_bma;
  job_data->ob_data.t[iteration] = dyn_ob->PE_ta;
  job_data->ob_data.um[iteration] = dyn_ob->PE_uma;

  if (iteration % 10000 == 0) {
    printf("PE calculation progress: %lld / %lld, %g%%                \r", iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }
  if (iteration == job_data->ob_data.len - 1) printf("Finished generating PE data for seed %f, process took %g seconds                \n", RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
}

void store_onebound_PEs_and_forces_callback(void* dyn, State s, void* job_msg, data_union *job_data, long long iteration) {
  assert(s == NEARBOUND);
  assert(iteration < job_data->ob_data.len);

  double* utility_ptr = (double*) job_msg;

  double start_time = *(utility_ptr[0]);
  long long max_iteration = *((long long*) utility_ptr[1]);
  
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
  job_data->ob_data.bb[iteration] = dyn_ob->PE_bba;
  job_data->ob_data.bm[iteration] = dyn_ob->PE_bma;
  job_data->ob_data.t[iteration] = dyn_ob->PE_ta;
  job_data->ob_data.um[iteration] = dyn_ob->PE_uma;

  utility_ptr[2][iteration] = f.bbx;
  utility_ptr[3][iteration] = f.bby;
  utility_ptr[4][iteration] = f.bmx;
  utility_ptr[5][iteration] = f.bmy;
  utility_ptr[6][iteration] = f.tx;
  utility_ptr[7][iteration] = f.ty;
  utility_ptr[8][iteration] = f.umx;
  utility_ptr[9][iteration] = f.umy;
  utility_ptr[10][iteration] = f.ubx;
  utility_ptr[11][iteration] = f.uby;

  if (iteration % 10000 == 0) {
    printf("PE calculation progress: %lld / %lld, %g%%                \r", iteration, max_iteration, ((double) iteration) / max_iteration * 100);
    fflush(NULL);
  }
  if (iteration == job_data->ob_data.len - 1) printf("Finished generating PE data for seed %f, process took %g seconds                \n", RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
}

void print_data_to_file(double* data1, double* data2, int iterations, const char* legend, const char* fname) {
  char* buf = (char*)
    malloc((iterations * (12+5+12+1) + 1 + 1) * sizeof(char) + strlen(legend)); // 5 spaced, each double 12 chars, newline, legend
  int buf_offset = 0;
  sprintf(buf, "%s\n", legend);
  buf_offset += strlen(legend) + 1;

  for (int i = 0; i < iterations; i++) {
    assert(data2[i] == data2[i]);
    sprintf(&buf[buf_offset], "%+.5e     %+.5e\n", data1[i], data2[i]);
    buf_offset += 30;
  }

  FILE* data_file = fopen(fname, "w");
  fputs(buf, data_file);
  fclose(data_file);
  free(buf);
}

void get_onebound_PE_correlation_function(generic_data* tau_data, onebound_data* corr_data, long long d_tau_iter, long long iterations, long long max_tau_iter, const int* seeds, int seed_len) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  long long num_corr_datapoints = max_tau_iter / d_tau_iter;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double test_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  onebound_data data;
  data.bb = (double*) malloc(iterations * sizeof(double));
  data.bm = (double*) malloc(iterations * sizeof(double));
  data.t  = (double*) malloc(iterations * sizeof(double));
  data.um = (double*) malloc(iterations * sizeof(double));
  data.len = iterations;

  data_union data_holder;
  data_holder.ob_data = data;
  double job_msg[2];
  job_msg[0] = (double) iterations;
  job_msg[1] = (double) clock();

  for (int i=0; i<num_corr_datapoints; i++) {
    corr_data->bb[i] = 0;
    corr_data->bm[i] = 0;
    corr_data->t[i]  = 0;
    corr_data->um[i] = 0;
  }
  
  for (int s = 0; s < seed_len; s++) {
    RAND_INIT_SEED = seeds[s];
    simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, test_position, store_onebound_PEs_callback, (void*) job_msg, &data_holder);

    double* PE_bbas = data.bb;
    double* PE_bmas = data.bm;
    double* PE_tas =  data.t;
    double* PE_umas = data.um;

    double PE_bba_ave = get_average(PE_bbas, iterations);
    double PE_bba_var = get_variance(PE_bbas, iterations);
    double PE_bma_ave = get_average(PE_bmas, iterations);
    double PE_bma_var = get_variance(PE_bmas, iterations);
    double PE_ta_ave = get_average(PE_tas, iterations);
    double PE_ta_var = get_variance(PE_tas, iterations);
    double PE_uma_ave = get_average(PE_umas, iterations);
    double PE_uma_var = get_variance(PE_umas, iterations);

    double start_time = (double) clock();

    for (long long tau_iter = 0; tau_iter < max_tau_iter; tau_iter += d_tau_iter) {
      double bba_correlation_sum = 0;
      double bma_correlation_sum = 0;
      double ta_correlation_sum = 0;
      double uma_correlation_sum = 0;
    
      for (long long i = 0; i < iterations - tau_iter; i++) {
	bba_correlation_sum += (PE_bbas[i] - PE_bba_ave) * (PE_bbas[i+tau_iter] - PE_bba_ave);
	bma_correlation_sum += (PE_bmas[i] - PE_bma_ave) * (PE_bmas[i+tau_iter] - PE_bma_ave);
	ta_correlation_sum += (PE_tas[i] - PE_ta_ave) * (PE_tas[i+tau_iter] - PE_ta_ave);
	uma_correlation_sum += (PE_umas[i] - PE_uma_ave) * (PE_umas[i+tau_iter] - PE_uma_ave);
      }

      long long num_iters = iterations - tau_iter;

      double bba_correlation = bba_correlation_sum / PE_bba_var / num_iters;
      double bma_correlation = bma_correlation_sum / PE_bma_var / num_iters;
      double ta_correlation  = ta_correlation_sum / PE_ta_var / num_iters;
      double uma_correlation = uma_correlation_sum / PE_uma_var / num_iters;

      long long iteration = tau_iter / d_tau_iter;
      assert(iteration < corr_data->len);

      if (tau_data != NULL) ((double*) tau_data->data)[iteration] = tau_iter*dt;
      corr_data->bb[iteration] += bba_correlation / seed_len; // average over seeds
      corr_data->bm[iteration] += bma_correlation / seed_len;
      corr_data->t[iteration]  += ta_correlation / seed_len;
      corr_data->um[iteration] += uma_correlation / seed_len;

      if (iteration % 1 == 0) {
	printf("correlation function progress: %g%%                \r", ((double) iteration)  / num_corr_datapoints * 100);
	fflush(NULL);
      }
      if (iteration == corr_data->len - 1) printf("Finished correlation calculations for seed %f, process took %g seconds                \n",
						  RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
  free(data.bb); free(data.bm); free(data.t); free(data.um);
}

void get_onebound_equipartition_ratio_per_runtime(generic_data* runtime_data, onebound_data* eq_data, long long d_runtime_iter, long long min_runtime_iter, long long max_runtime_iter, const int* seeds, int seed_len) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  long long num_eq_datapoints;
  if (min_runtime_iter == max_runtime_iter)
    num_eq_datapoints = 1;
  else
     num_eq_datapoints = (max_runtime_iter - min_runtime_iter) / d_runtime_iter;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  onebound_data data;
  data.bb = (double*) malloc(max_runtime_iter * sizeof(double));
  data.bm = (double*) malloc(max_runtime_iter * sizeof(double));
  data.t  = (double*) malloc(max_runtime_iter * sizeof(double));
  data.um = (double*) malloc(max_runtime_iter * sizeof(double));
  data.len = max_runtime_iter;

  data_union data_holder;
  data_holder.ob_data = data;

  for (int i=0; i<num_eq_datapoints; i++) {
    eq_data->bb[i] = 0;
    eq_data->bm[i] = 0;
    eq_data->t[i]  = 0;
    eq_data->um[i] = 0;
  }
  
  for (int s = 0; s < seed_len; s++) {
    RAND_INIT_SEED = seeds[s];

    double job_msg[2];
    job_msg[0] = (double) max_runtime_iter;
    job_msg[1] = (double) clock();
    
    simulate(max_runtime_iter*dt, RAND_INIT_SEED, NEARBOUND, init_position, store_onebound_PEs_callback, (void*) job_msg, &data_holder);

    double* bba_PE_data = data_holder.ob_data.bb;
    double* bma_PE_data = data_holder.ob_data.bm;
    double* ta_PE_data =  data_holder.ob_data.t; 
    double* uma_PE_data = data_holder.ob_data.um;

    double start_time = (double) clock();

    for (long long runtime_iter = min_runtime_iter; runtime_iter < max_runtime_iter; runtime_iter += d_runtime_iter) {
      long long iteration = (runtime_iter - min_runtime_iter) / d_runtime_iter;
      assert(iteration < eq_data->len);

      if (runtime_data != NULL) ((double*) runtime_data->data)[iteration] = runtime_iter * dt;
      eq_data->bb[iteration] += bba_PE_data[runtime_iter] / (0.5*kb*T) / seed_len; // seed average
      eq_data->bm[iteration] += bma_PE_data[runtime_iter] / (0.5*kb*T) / seed_len;
      eq_data->t[iteration]  += ta_PE_data[runtime_iter] / (0.5*kb*T)  / seed_len;
      eq_data->um[iteration] += uma_PE_data[runtime_iter] / (0.5*kb*T) / seed_len;

      if (runtime_iter/d_runtime_iter % 1 == 0) {
	printf("equipartition ratio progress: %g%%                \r", ((double) iteration) / num_eq_datapoints * 100);
	fflush(NULL);
      }
      if (iteration == eq_data->len - 1) printf("Finished equipartition calculations for seed %f, process took %g seconds                \n",
						RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
  free(data.bb); free(data.bm); free(data.t); free(data.um);
}

void get_onebound_equipartition_ratio_average_per_runtime(generic_data* runtime_data, onebound_data* eq_data, long long d_runtime_iter, long long min_runtime_iter, long long max_runtime_iter, const int* seeds, int seed_len) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  long long num_eq_datapoints;
  if (min_runtime_iter == max_runtime_iter)
    num_eq_datapoints = 1;
  else
     num_eq_datapoints = (max_runtime_iter - min_runtime_iter) / d_runtime_iter;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  onebound_data data;
  data.bb = (double*) malloc(max_runtime_iter * sizeof(double));
  data.bm = (double*) malloc(max_runtime_iter * sizeof(double));
  data.t  = (double*) malloc(max_runtime_iter * sizeof(double));
  data.um = (double*) malloc(max_runtime_iter * sizeof(double));
  data.len = max_runtime_iter;

  data_union data_holder;
  data_holder.ob_data = data;

  for (int i=0; i<num_eq_datapoints; i++) {
    eq_data->bb[i] = 0;
    eq_data->bm[i] = 0;
    eq_data->t[i]  = 0;
    eq_data->um[i] = 0;
  }
  
  for (int s = 0; s < seed_len; s++) {
    RAND_INIT_SEED = seeds[s];

    double job_msg[2];
    job_msg[0] = (double) max_runtime_iter;
    job_msg[1] = (double) clock();

    simulate(max_runtime_iter*dt, RAND_INIT_SEED, NEARBOUND, init_position, store_onebound_PEs_callback, (void*) job_msg, &data_holder);

    double* bba_PE_data = data_holder.ob_data.bb;
    double* bma_PE_data = data_holder.ob_data.bm;
    double* ta_PE_data =  data_holder.ob_data.t; 
    double* uma_PE_data = data_holder.ob_data.um;

    double start_time = (double) clock();

    for (long long runtime_iter = min_runtime_iter; runtime_iter < max_runtime_iter; runtime_iter += d_runtime_iter) {
      double bba_eq_ratio = get_average(bba_PE_data, runtime_iter) / (0.5*kb*T);
      double bma_eq_ratio = get_average(bma_PE_data, runtime_iter) / (0.5*kb*T);
      double  ta_eq_ratio = get_average( ta_PE_data, runtime_iter) / (0.5*kb*T);
      double uma_eq_ratio = get_average(uma_PE_data, runtime_iter) / (0.5*kb*T);

      long long iteration = (runtime_iter - min_runtime_iter) / d_runtime_iter;
      assert(iteration < eq_data->len);

      if (runtime_data != NULL) ((double*) runtime_data->data)[iteration] = runtime_iter * dt;
      eq_data->bb[iteration] += bba_eq_ratio / seed_len;
      eq_data->bm[iteration] += bma_eq_ratio / seed_len;
      eq_data->t[iteration]  += ta_eq_ratio  / seed_len;
      eq_data->um[iteration] += uma_eq_ratio / seed_len;

      if (runtime_iter/d_runtime_iter % 1 == 0) {
	printf("equipartition ratio progress: %g%%                \r", ((double) iteration) / num_eq_datapoints * 100);
	fflush(NULL);
      }
      if (iteration == eq_data->len - 1) printf("Finished equipartition calculations for seed %f, process took %g seconds                \n",
						RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
  free(data.bb); free(data.bm); free(data.t); free(data.um);
}

void get_onebound_equipartition_ratio_and_average_force(generic_data* force_data, onebound_data* eq_data, long long runtime_iter, const int* seeds, int seed_len) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  onebound_data data;
  data.bb = (double*) malloc(runtime_iter * sizeof(double));
  data.bm = (double*) malloc(runtime_iter * sizeof(double));
  data.t  = (double*) malloc(runtime_iter * sizeof(double));
  data.um = (double*) malloc(runtime_iter * sizeof(double));
  data.len = max_runtime_iter;

  data_union data_holder;
  data_holder.ob_data = data;

  eq_data->bb = 0;
  eq_data->bm = 0;
  eq_data->t  = 0;
  eq_data->um = 0;
  
  ((onebound_forces*) force_data->data)->bbx = 0;
  ((onebound_forces*) force_data->data)->bby = 0;
  ((onebound_forces*) force_data->data)->bmx = 0;
  ((onebound_forces*) force_data->data)->bmy = 0;
  ((onebound_forces*) force_data->data)->tx  = 0;
  ((onebound_forces*) force_data->data)->ty  = 0;
  ((onebound_forces*) force_data->data)->umx = 0;
  ((onebound_forces*) force_data->data)->umy = 0;
  ((onebound_forces*) force_data->data)->ubx = 0;
  ((onebound_forces*) force_data->data)->uby = 0;

  double job_msg[12];
  job_msg[0] = (double*) &max_runtime_iter;
  job_msg[2] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[3] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[4] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[5] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[6] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[7] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[8] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[9] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[10] = (double*) malloc(runtime_iter * sizeof(double));
  job_msg[11] = (double*) malloc(runtime_iter * sizeof(double));
  
  for (int s = 0; s < seed_len; s++) {
    RAND_INIT_SEED = seeds[s];

    double current_time = clock();
    job_msg[1] = (double*) &current_time;
    
    simulate(max_runtime_iter*dt, RAND_INIT_SEED, NEARBOUND, init_position,
	store_onebound_PEs_and_forces_callback, (void*) job_msg, &data_holder);
    
    double bba_eq_ratio = get_average(data.bb, runtime_iter) / (0.5*kb*T);
    double bma_eq_ratio = get_average(data.bm, runtime_iter) / (0.5*kb*T);
    double  ta_eq_ratio = get_average(data.t,  runtime_iter) / (0.5*kb*T);
    double uma_eq_ratio = get_average(data.um, runtime_iter) / (0.5*kb*T);

    eq_data->bb += bba_eq_ratio / seed_len;
    eq_data->bm += bma_eq_ratio / seed_len;
    eq_data->t  += ta_eq_ratio  / seed_len;
    eq_data->um += uma_eq_ratio / seed_len;

    double f_bbx_variance = get_variance(job_msg[2], runtime_iter);
    double f_bby_variance = get_variance(job_msg[3], runtime_iter);
    double f_bmx_variance = get_variance(job_msg[4], runtime_iter);
    double f_bmy_variance = get_variance(job_msg[5], runtime_iter);
    double f_tx_variance  = get_variance(job_msg[6], runtime_iter);
    double f_ty_variance  = get_variance(job_msg[7], runtime_iter);
    double f_umx_variance = get_variance(job_msg[8], runtime_iter);
    double f_umy_variance = get_variance(job_msg[9], runtime_iter);
    double f_ubx_variance = get_variance(job_msg[10], runtime_iter);
    double f_uby_variance = get_variance(job_msg[11], runtime_iter);

    ((onebound_forces*) force_data->data)->bbx += f_bbx_variance / seed_len;
    ((onebound_forces*) force_data->data)->bby += f_bby_variance / seed_len;
    ((onebound_forces*) force_data->data)->bmx += f_bmx_variance / seed_len;
    ((onebound_forces*) force_data->data)->bmy += f_bmy_variance / seed_len;
    ((onebound_forces*) force_data->data)->tx  += f_tx_variance / seed_len;
    ((onebound_forces*) force_data->data)->ty  += f_ty_variance / seed_len;
    ((onebound_forces*) force_data->data)->umx += f_umx_variance / seed_len;
    ((onebound_forces*) force_data->data)->umy += f_umy_variance / seed_len;
    ((onebound_forces*) force_data->data)->ubx += f_ubx_variance / seed_len;
    ((onebound_forces*) force_data->data)->uby += f_uby_variance / seed_len;
  }

  free(data.bb); free(data.bm); free(data.t); free(data.um);
  free(job_msg[2]); free(job_msg[3]); free(job_msg[4]);
  free(job_msg[5]); free(job_msg[6]); free(job_msg[7]);
  free(job_msg[8]); free(job_msg[9]); free(job_msg[10]);
  free(job_msg[11]);
}
