#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "../dynein_struct.h"
#include "simulation_defaults.h"

/** Library for simulation code **/

void store_onebound_PEs_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  if (iteration % data_generation_skip_iterations != 0) return;
  else {
    assert(s == NEARBOUND);
    assert(iteration < job_data->ob_data.len);
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];

    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
    job_data->ob_data.bb[iteration] = dyn_ob->PE_bba;
    job_data->ob_data.bm[iteration] = dyn_ob->PE_bma;
    job_data->ob_data.t[iteration] = dyn_ob->PE_ta;
    job_data->ob_data.um[iteration] = dyn_ob->PE_uma;

    if (iteration % 10000 == 0) {
      printf("PE calculation progress for %s: %lld / %lld, %g%%                \r",
	     run_msg, iteration, max_iteration, ((double) iteration) / max_iteration * 100);
      fflush(NULL);
    }
    if (iteration == job_data->ob_data.len - 1)
      printf("Finished generating PE data for %s, which took %g seconds                \n",
	     run_msg, ((double) clock() - start_time) / CLOCKS_PER_SEC);
  }
}

void store_onebound_PEs_and_forces_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  if (iteration % data_generation_skip_iterations != 0) return;
  else {
    assert(s == NEARBOUND);
    int iter = iteration / data_generation_skip_iterations;

    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];

    eq_forces_callback_data_ptr data = *((eq_forces_callback_data_ptr*) job_data->g_data.data);
    assert(iter < job_data->g_data.len);

    data.bb[iter] = dyn_ob->PE_bba;
    data.bm[iter] = dyn_ob->PE_bma;
    data.t[iter]  = dyn_ob->PE_ta;
    data.um[iter] = dyn_ob->PE_uma;

    onebound_forces f = dyn_ob->get_internal();

    if (data.f_bbx != NULL) {
      data.f_bbx[iter] = f.bbx;
      data.f_bby[iter] = f.bby;
      data.f_bmx[iter] = f.bmx;
      data.f_bmy[iter] = f.bmy;
      data.f_tx [iter] = f.tx;
      data.f_ty [iter] = f.ty;
      data.f_umx[iter] = f.umx;
      data.f_umy[iter] = f.umy;
      data.f_ubx[iter] = f.ubx;
      data.f_uby[iter] = f.uby;
    }

    if (iter % 10 == 0) {
      printf("PE calculation progress for %s: %d / %lld, %g%%                \r",
	     run_msg, iter, max_iteration, ((double) iter) / max_iteration * 100);
      fflush(NULL);
    }
    if (iter == job_data->g_data.len - 1)
      printf("Finished generating PE data for %s, which took %g seconds                \n",
	     run_msg, ((double) clock() - start_time) / CLOCKS_PER_SEC);
  }
}


void get_onebound_PE_correlation_function(generic_data* tau_data, onebound_data* corr_data, long long num_correlations, long long iterations, long long max_tau_iter, const int* seeds, int seed_len, char* run_msg_base) {
  long long d_tau_iter = max_tau_iter / num_correlations;

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

  void* job_msg[3];
  job_msg[0] = (double*) &iterations;
  double current_time = clock();
  job_msg[1] = &current_time;

  char run_msg[512];
  strcpy(run_msg, run_msg_base);
  char seedbuf[50];
  sprintf(seedbuf, "seed = %d)", (int) RAND_INIT_SEED);
  strcat(run_msg, seedbuf);
  job_msg[2] = run_msg;

  for (int i=0; i<num_correlations; i++) {
    corr_data->bb[i] = 0;
    corr_data->bm[i] = 0;
    corr_data->t[i]  = 0;
    corr_data->um[i] = 0;
  }

  for (int s = 0; s < seed_len; s++) {
    RAND_INIT_SEED = seeds[s];
    simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, test_position,
	     store_onebound_PEs_callback, job_msg, &data_holder);

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
	printf("correlation function progress: %g%%                \r", ((double) iteration)  / num_correlations * 100);
	fflush(NULL);
      }
      if (iteration == corr_data->len - 1) printf("Finished correlation calculations for seed %f, process took %g seconds                \n",
						  RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
  free(data.bb); free(data.bm); free(data.t); free(data.um);
}

void get_onebound_equipartition_ratio_per_runtime(generic_data* runtime_data, onebound_data* eq_data, long long num_eq_datapoints, long long min_runtime_iter, long long max_runtime_iter, const int* seeds, int seed_len, char* run_msg_base) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;

  long long d_runtime_iter = (max_runtime_iter - min_runtime_iter) / num_eq_datapoints;

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

    void* job_msg[3];
    job_msg[0] = (double*) &max_runtime_iter;
    double current_time = clock();
    job_msg[1] = &current_time;
    
    char run_msg[512];
    strcpy(run_msg, run_msg_base);
    char seedbuf[50];
    sprintf(seedbuf, "seed = %d)", (int) RAND_INIT_SEED);
    strcat(run_msg, seedbuf);
    job_msg[2] = run_msg;

    simulate(max_runtime_iter*dt, RAND_INIT_SEED, NEARBOUND, init_position, store_onebound_PEs_callback, job_msg, &data_holder);

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

void get_onebound_equipartition_ratio_average_per_runtime(generic_data* runtime_data, onebound_data* eq_data, long long num_eq_datapoints, long long min_runtime_iter, long long max_runtime_iter, const int* seeds, int seed_len, char* run_msg_base) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;

  long long d_runtime_iter = (max_runtime_iter - min_runtime_iter) / num_eq_datapoints;

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

    void* job_msg[3];
    job_msg[0] = (double*) &max_runtime_iter;
    double current_time = clock();
    job_msg[1] = &current_time;
    
    char run_msg[512];
    strcpy(run_msg, run_msg_base);
    char seedbuf[50];
    sprintf(seedbuf, "seed = %d)", (int) RAND_INIT_SEED);
    strcat(run_msg, seedbuf);
    job_msg[2] = run_msg;

    simulate(max_runtime_iter*dt, RAND_INIT_SEED, NEARBOUND, init_position, store_onebound_PEs_callback, job_msg, &data_holder);

    double* bba_PE_data = data_holder.ob_data.bb;
    double* bma_PE_data = data_holder.ob_data.bm;
    double* ta_PE_data =  data_holder.ob_data.t; 
    double* uma_PE_data = data_holder.ob_data.um;

    double start_time = (double) clock();

    for (long long i = 0; i < num_eq_datapoints; i++) {
      long long runtime_iter = min_runtime_iter + i*d_runtime_iter;
      double bba_eq_ratio = get_average(bba_PE_data, runtime_iter) / (0.5*kb*T);
      double bma_eq_ratio = get_average(bma_PE_data, runtime_iter) / (0.5*kb*T);
      double  ta_eq_ratio = get_average( ta_PE_data, runtime_iter) / (0.5*kb*T);
      double uma_eq_ratio = get_average(uma_PE_data, runtime_iter) / (0.5*kb*T);

      assert(i < eq_data->len);

      if (runtime_data != NULL) ((double*) runtime_data->data)[i] = runtime_iter * dt;
      eq_data->bb[i] += bba_eq_ratio / seed_len;
      eq_data->bm[i] += bma_eq_ratio / seed_len;
      eq_data->t[i]  += ta_eq_ratio  / seed_len;
      eq_data->um[i] += uma_eq_ratio / seed_len;

      if (runtime_iter/d_runtime_iter % 1 == 0) {
	printf("equipartition ratio progress: %g%%                \r",
	     ((double) i) / num_eq_datapoints * 100);
	fflush(NULL);
      }
      if (i == eq_data->len - 1) printf("Finished equipartition calculations for seed %f,"
	"process took %g seconds                \n",
        RAND_INIT_SEED, ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
  free(data.bb); free(data.bm); free(data.t); free(data.um);
}

void get_onebound_equipartition_ratio(onebound_data* eq_data, generic_data* force_data, long long iters, const int* seeds, int seed_len, char* run_msg_base) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba,
			    eq.bma + eq.bba - M_PI,
			    eq.ta + eq.bma + eq.bba - M_PI,
			    eq.ta + eq.bma + eq.bba - eq.uma,
			    0, 0};

  eq_forces_callback_data_ptr data_ptr;
  generic_data data_struct;
  data_struct.data = &data_ptr;
  data_struct.len = iters;

  data_union data_holder;
  data_holder.g_data = data_struct;

  data_ptr.bb = (double*) malloc(iters * sizeof(double));
  data_ptr.bm = (double*) malloc(iters * sizeof(double));
  data_ptr.t  = (double*) malloc(iters * sizeof(double));
  data_ptr.um = (double*) malloc(iters * sizeof(double));

  eq_data->bb[0] = 0;    eq_data->t[0]  = 0;
  eq_data->bm[0] = 0;    eq_data->um[0] = 0;

  data_ptr.f_bbx = (double*) malloc(iters * sizeof(double));
  data_ptr.f_bby = (double*) malloc(iters * sizeof(double));
  data_ptr.f_bmx = (double*) malloc(iters * sizeof(double));
  data_ptr.f_bmy = (double*) malloc(iters * sizeof(double));
  data_ptr.f_tx  = (double*) malloc(iters * sizeof(double));
  data_ptr.f_ty  = (double*) malloc(iters * sizeof(double));
  data_ptr.f_umx = (double*) malloc(iters * sizeof(double));
  data_ptr.f_umy = (double*) malloc(iters * sizeof(double));
  data_ptr.f_ubx = (double*) malloc(iters * sizeof(double));
  data_ptr.f_uby = (double*) malloc(iters * sizeof(double));

  double* f_bbx_var = &((onebound_forces*) force_data->data)->bbx;
  double* f_bby_var = &((onebound_forces*) force_data->data)->bby;
  double* f_bmx_var = &((onebound_forces*) force_data->data)->bmx;
  double* f_bmy_var = &((onebound_forces*) force_data->data)->bmy;
  double* f_tx_var  = &((onebound_forces*) force_data->data)->tx; 
  double* f_ty_var  = &((onebound_forces*) force_data->data)->ty; 
  double* f_umx_var = &((onebound_forces*) force_data->data)->umx;
  double* f_umy_var = &((onebound_forces*) force_data->data)->umy;
  double* f_ubx_var = &((onebound_forces*) force_data->data)->ubx;
  double* f_uby_var = &((onebound_forces*) force_data->data)->uby;

  *f_bbx_var = 0;   *f_bby_var = 0;   *f_bmx_var = 0;   *f_bmy_var = 0;
  *f_tx_var = 0;    *f_ty_var = 0;    *f_umx_var = 0;   *f_umy_var = 0;
  *f_ubx_var = 0;   *f_uby_var = 0;

  void* job_msg[3];
  job_msg[0] = (double*) &iters;

  char run_msg[512];

  for (int s = 0; s < seed_len; s++) {
    RAND_INIT_SEED = seeds[s];

    double current_time = clock();
    job_msg[1] = &current_time;

    strcpy(run_msg, run_msg_base);
    char seedbuf[50];
    sprintf(seedbuf, "seed = %d)", (int) RAND_INIT_SEED);
    strcat(run_msg, seedbuf);
    job_msg[2] = run_msg;

    simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, init_position,
	store_onebound_PEs_and_forces_callback, job_msg, &data_holder);

    eq_data->bb[0] += get_average(data_ptr.bb, iters) / (0.5*kb*T) / seed_len;
    eq_data->bm[0] += get_average(data_ptr.bm, iters) / (0.5*kb*T) / seed_len;
    eq_data->t [0] += get_average(data_ptr.t,  iters) / (0.5*kb*T) / seed_len;
    eq_data->um[0] += get_average(data_ptr.um, iters) / (0.5*kb*T) / seed_len;

    *f_bbx_var += get_variance(data_ptr.f_bbx, iters) / seed_len;
    *f_bby_var += get_variance(data_ptr.f_bby, iters) / seed_len;
    *f_bmx_var += get_variance(data_ptr.f_bmx, iters) / seed_len;
    *f_bmy_var += get_variance(data_ptr.f_bmy, iters) / seed_len;
    *f_tx_var  += get_variance(data_ptr.f_tx , iters) / seed_len;
    *f_ty_var  += get_variance(data_ptr.f_ty , iters) / seed_len;
    *f_umx_var += get_variance(data_ptr.f_umx, iters) / seed_len;
    *f_umy_var += get_variance(data_ptr.f_umy, iters) / seed_len;
    *f_ubx_var += get_variance(data_ptr.f_ubx, iters) / seed_len;
    *f_uby_var += get_variance(data_ptr.f_uby, iters) / seed_len;
  }

  free(data_ptr.bb);  free(data_ptr.bm);  free(data_ptr.t);  free(data_ptr.um);
  free(data_ptr.f_bbx);  free(data_ptr.f_bby);  free(data_ptr.f_bmx);
  free(data_ptr.f_bmy);  free(data_ptr.f_tx );  free(data_ptr.f_ty );
  free(data_ptr.f_umx);  free(data_ptr.f_umy);  free(data_ptr.f_ubx);
  free(data_ptr.f_uby);
}
