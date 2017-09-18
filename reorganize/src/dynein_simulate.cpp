#include "dynein_struct.h"
#include "simulation_defaults.h"
#include <stdlib.h> // for exit
#include <string.h>

static const bool debug_stepping = false;

double variable_ts_checkpoint_interval = 1e-6;
char* variable_ts_stepping_print_buffer;
int variable_ts_stepping_print_buffer_index;
FILE* variable_ts_stepping_data_file;

Dynein_onebound* variable_ts_checkpoint_onebound;
Dynein_bothbound* variable_ts_checkpoint_bothbound;
double variable_ts_checkpoint_time;

int variable_ts_rewinding_state = false;
double variable_ts_base_dt;

void simulate(double runtime, double rand_seed, State init_state, double* init_position,
	      void (*job)(void* dyn, State s, void *job_msg, data_union* job_data,
	      long long iteration), void *job_msg, data_union* job_data) {

  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT);      // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  
  MTRand* rand = new MTRand(rand_seed);

  Dynein_onebound *dyn_ob;
  Dynein_bothbound *dyn_bb;
  
  if (init_state == BOTHBOUND) {
    dyn_ob = NULL;
    dyn_bb = new Dynein_bothbound(
			      init_position[0],      // nma_init
				  init_position[1],      // fma_init
				  init_position[2],      // nbx_init
				  init_position[3],      // nby_init
				  init_position[4],      // L
				  NULL,                  // internal forces
				  NULL,                  // brownian forces
				  NULL,                  // equilibrium angles
				  rand);                 // MTRand
  } else {
    dyn_ob = new Dynein_onebound(
				 init_position[0],    // bba_init
				 init_position[1],    // bma_init
				 init_position[2],    // uma_init
				 init_position[3],    // uba_init
				 init_position[4],    // bbx_init
				 init_position[5],    // bby_init
				 init_state,          // Initial state
				 NULL,                // Optional custom internal forces
				 NULL,                // Optional custom brownian forces
				 NULL,                // Optional custom equilibrium angles
				 rand);
    dyn_bb = NULL;
  }

  double t = 0;
  long long iter = 0;
  State current_state = init_state;

  double rebinding_immune_until = 0; // to prevent immediate rebinding in BB->OB transitions

  bool run_indefinite;
  if (runtime == 0) {
	run_indefinite = true;
    printf("Running indefinitely.\n");
  } else {
	run_indefinite = false;
	printf("Running for %g s\n", runtime);
  }

  double near_unbinding_prob_printing_average = 0;
  int unbinding_prob_printing_n = 0;

  while( t < runtime or run_indefinite) {
    if (current_state == NEARBOUND or current_state == FARBOUND)
      while (t < runtime or run_indefinite) { // loop as long as it is onebound
        if (am_debugging_time) printf("\n==== t = %8g/%8g ====\n", t, runtime);

        double unbinding_prob = dyn_ob->get_unbinding_rate()*dt;
        double binding_prob = dyn_ob->get_binding_rate()*dt;
	if (am_debugging_rates and binding_prob != 0 and rand->rand() < 1e-3) {
	  printf("binding probability: %g, uby %g at time %g s\n", binding_prob, dyn_ob->get_uby(), t);
	}
	if (rand->rand() < unbinding_prob) { // unbind, switch to unbound
	  delete dyn_ob;
	  dyn_ob = NULL;
	  current_state = UNBOUND;
	  break;
	}
	else if (rand->rand() < binding_prob and t > rebinding_immune_until) { // switch to bothbound
	  dyn_bb = new Dynein_bothbound(dyn_ob, rand);
	  if (am_debugging_state_transitions) printf("Transitioning from onebound to bothbound\n");
	  if (am_debugging_state_transitions) printf("just bound b/c binding probability was: %g, boltzmann factor: %g\n",
                                                     binding_prob, exp(-(dyn_bb->get_PE()-dyn_ob->get_PE())/kb/T));
	  delete dyn_ob;
	  dyn_ob = NULL;
	  current_state = BOTHBOUND;
	  job(dyn_bb, current_state, job_msg, job_data, iter);
	  t += dt;
	  iter++;
	  break;
	}
	else { // move like normal
	  job(dyn_ob, current_state, job_msg, job_data, iter);
	  t += dt;
	  iter++;

	  // potentially faster to compute velocity here, instead of down there?

	  double temp_bba = dyn_ob->get_bba() + dyn_ob->get_d_bba() * dt;
	  double temp_bma = dyn_ob->get_bma() + dyn_ob->get_d_bma() * dt;
	  double temp_uma = dyn_ob->get_uma() + dyn_ob->get_d_uma() * dt;
	  double temp_uba = dyn_ob->get_uba() + dyn_ob->get_d_uba() * dt;

	  dyn_ob->set_bba(temp_bba);
	  dyn_ob->set_bma(temp_bma);
	  dyn_ob->set_uma(temp_uma);
	  dyn_ob->set_uba(temp_uba);

	  dyn_ob->update_velocities();

	  if (crash_on_nan and (isnan(temp_bba) or isnan(temp_bma) or isnan(temp_uma) or isnan(temp_uba))) {
	    printf("Onebound velocity calculation generated a NaN, exiting.\n");
	    exit(1);
	  }
	}
      }

    if (current_state == BOTHBOUND) {
      while (t < runtime or run_indefinite) { // loop as long as it is bothbound

        if (am_debugging_time) printf("\n==== t = %8g/%8g ====\n", t, runtime);
        double near_unbinding_prob = dyn_bb->get_near_unbinding_rate()*dt;
        double far_unbinding_prob = dyn_bb->get_far_unbinding_rate()*dt;
	double roll = rand->rand();
	while (roll == 0) roll = rand->rand();
	if (am_debugging_rates and roll < 1e-8) {
	  near_unbinding_prob_printing_average = (near_unbinding_prob_printing_average*unbinding_prob_printing_n + near_unbinding_prob);
	  near_unbinding_prob_printing_average /= (unbinding_prob_printing_n + 1);
	  unbinding_prob_printing_n++;

	  printf("BB near unbinding probability: %g\n", near_unbinding_prob_printing_average);
	}
	bool unbind_near = roll < near_unbinding_prob;
	bool unbind_far = roll < far_unbinding_prob;
	if (unbind_near && unbind_far) {
	  if (debug_stepping) printf("both MTBDs want to fall off!\n");
	  if (iter % 2 == 0) unbind_far = false;
	  else unbind_near = false;
	}
	if (unbind_near) { // switch to farbound
	  dyn_ob = new Dynein_onebound(dyn_bb, rand, FARBOUND);
	  if (am_debugging_state_transitions) printf("Transitioning from bothbound to farbound\n");
	  if (am_debugging_state_transitions) printf("just unbound b/c unbinding probability was: %g, roll was: %g, boltzmann factor: %g\n",
                                                     near_unbinding_prob, roll, exp(-(dyn_ob->get_PE()-dyn_bb->get_PE())/kb/T));
	  delete dyn_bb;
	  dyn_bb = NULL;
	  current_state = FARBOUND;
	  job(dyn_ob, current_state, job_msg, job_data, iter);
	  t += dt;
	  iter++;
	  rebinding_immune_until = t + REBINDING_IMMUNITY_TIME;
	  break;
	}
	else if (unbind_far) { // switch to nearbound
	  dyn_ob = new Dynein_onebound(dyn_bb, rand, NEARBOUND);
	  if (am_debugging_state_transitions) printf("Transitioning from bothbound to nearbound\n");
	  if (am_debugging_state_transitions) printf("just unbound b/c unbinding probability was: %g, roll as: %g, boltzmann factor: %g\n",
                                                     far_unbinding_prob, roll, exp(-(dyn_ob->get_PE()-dyn_bb->get_PE())/kb/T));
	  delete dyn_bb;
	  dyn_bb = NULL;
	  current_state = NEARBOUND;
	  job(dyn_ob, current_state, job_msg, job_data, iter);
	  t += dt;
	  iter++;
	  rebinding_immune_until = t + REBINDING_IMMUNITY_TIME;
	  break;
	}
	else { // move like normal
	  job(dyn_bb, BOTHBOUND, job_msg, job_data, iter);
	  t += dt;
	  iter++;

	  double temp_nma = dyn_bb->get_nma() + dyn_bb->get_d_nma()*dt;
	  double temp_fma = dyn_bb->get_fma() + dyn_bb->get_d_fma()*dt;
	  dyn_bb->set_nma(temp_nma);
	  dyn_bb->set_fma(temp_fma);

	  dyn_bb->update_velocities();

	  if (crash_on_nan and (isnan(temp_nma) or isnan(temp_fma))) {
	    printf("Bothbound velocity calculation generated a NaN, exiting.\n");
	    exit(1);
	  }
	}
      }
    }
    if (current_state == UNBOUND) {
      job(NULL, UNBOUND, job_msg, job_data, iter);
      goto end_simulation;
    }
  }

 end_simulation:
  delete rand;
  if (dyn_bb == NULL) delete dyn_ob;
  else delete dyn_bb;
}
