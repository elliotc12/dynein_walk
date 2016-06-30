#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

void store_PE_bmas(void* dyn, State s, void* job_msg, void* job_data, int iteration) {
  assert(s == NEARBOUND);
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
  ((double*) job_data)[iteration] = dyn_ob->PE_bba;
}

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  
  T = 100;
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  
  int iterations = 1e6;
  double runtime = dt*iterations;
  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double test_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};
  void* data = malloc(iterations * sizeof(double));
  
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, test_position, store_PE_bmas, NULL, data);

  /** compute correlation function for PE_nma **/

  double* PE_nmas = (double*) data;
  
  double PE_nma_ave = get_average(PE_nmas, iterations);
  double PE_nma_var = get_variance(PE_nmas, iterations);

  int d_iter = 1e2;
  int max_tau_iter = floor(iterations * 0.5);

  for (int tau_iter=0; tau_iter < max_tau_iter; tau_iter += d_iter) {
    double correlation_sum = 0;
    for (int i = 0; i < iterations - tau_iter; i++) {
      correlation_sum += (PE_nmas[i] - PE_nma_ave) * (PE_nmas[i+tau_iter] - PE_nma_ave);
    }
    int num_iterations = iterations - tau_iter;
    double correlation = correlation_sum / PE_nma_var / num_iterations;
    printf("correlation at tau_iter = %d is: %g\n", tau_iter, correlation);
  }
}
