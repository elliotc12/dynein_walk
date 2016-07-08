#include "dynein_struct.h"

void simulate(double runtime, double rand_seed, State init_state, double* init_position,
	      void (*job)(void* dyn, State s, void* job_msg, data_union* job_data, long long iteration),
	      void* job_msg, data_union* job_data) {

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
  State current_state = init_state; 

  while( t < runtime ) {
    if (current_state == NEARBOUND or current_state == FARBOUND)
      while (t < runtime) { // loop as long as it is onebound
	if (rand->rand() < dyn_ob->get_unbinding_rate()*dt) { // unbind, switch to unbound
	  printf("unbinding at %.1f%%!\n", t/runtime*100);
	  delete dyn_ob;
	  current_state = UNBOUND;
	  break;
	}
	else if (rand->rand() < dyn_ob->get_binding_rate()*dt) { // bind ub, switch to bothbound
	  printf("switch to bothbound at %.1f%%!\n", t/runtime*100);
	  dyn_bb = new Dynein_bothbound(dyn_ob, rand);
	  delete dyn_ob;
	  current_state = BOTHBOUND;
	  break;
	}
	else { // move like normal
	  double temp_bba = dyn_ob->get_bba() + dyn_ob->get_d_bba() * dt;
	  double temp_bma = dyn_ob->get_bma() + dyn_ob->get_d_bma() * dt;
	  double temp_uma = dyn_ob->get_uma() + dyn_ob->get_d_uma() * dt;
	  double temp_uba = dyn_ob->get_uba() + dyn_ob->get_d_uba() * dt;

	  dyn_ob->set_bba(temp_bba);
	  dyn_ob->set_bma(temp_bma);
	  dyn_ob->set_uma(temp_uma);
	  dyn_ob->set_uba(temp_uba);

	  long long iter = t/dt;
	  job(dyn_ob, current_state, job_msg, job_data, iter);
	  t += dt;
	  dyn_ob->update_velocities();
	}
      }

    if (current_state == BOTHBOUND) {
      while (t < runtime) { // loop as long as it is bothbound
	bool unbind_near = rand->rand() < dyn_bb->get_near_unbinding_rate();
	bool unbind_far = rand->rand() < dyn_bb->get_far_unbinding_rate();
	if (unbind_near && unbind_far) {
	  printf("THEY BOTH WANT TO FALL OFF TOGETHER!!!\n");
	  printf("WHAT SHOULD THEY DO????\n");
	}
	else if (unbind_near) { // switch to farbound
	  dyn_ob = new Dynein_onebound(dyn_bb, rand, FARBOUND);
	  delete dyn_bb;
	  printf("switch to onebound!\n");
	  break;
	}
	else if (unbind_far) { // switch to nearbound
	  dyn_ob = new Dynein_onebound(dyn_bb, rand, NEARBOUND);
	  delete dyn_bb;
	  printf("switch to onebound!\n");
	  break;
	}
	else { // move like normal
	  double temp_nma = dyn_bb->get_nma() + dyn_bb->get_d_nma()*dt;
	  double temp_fma = dyn_bb->get_fma() + dyn_bb->get_d_fma()*dt;
	  dyn_bb->set_nma(temp_nma);
	  dyn_bb->set_fma(temp_fma);

	  job(dyn_bb, BOTHBOUND, job_msg, job_data, (long long) t/dt);
	  t += dt;
	  dyn_bb->update_velocities();
	}
      }
    }
    if (current_state == UNBOUND) {
      while (t < runtime) {
	job(NULL, UNBOUND, job_msg, job_data, (long long) t/dt);
	t += dt;
      }
    }
  }
}
