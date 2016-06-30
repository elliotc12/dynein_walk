#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  
  int iterations = 1e7;
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

  int num_correlations = floor(max_tau_iter/d_iter);

  const char* correlation_legend = "PE_nma correlation function";
  char* correlation_buf = (char*)
    malloc((num_correlations * (12 + 5 + 12 + 1) + 1 + 1) * sizeof(char) + strlen(correlation_legend));
  //each line + legend + \n + \0
  int correlation_buf_offset = 0;
  sprintf(correlation_buf, "%s\n", correlation_legend);
  correlation_buf_offset += strlen(correlation_legend) + 1;

  for (int tau_iter=0; tau_iter < max_tau_iter; tau_iter += d_iter) {
    double correlation_sum = 0;
    for (int i = 0; i < iterations - tau_iter; i++) {
      correlation_sum += (PE_nmas[i] - PE_nma_ave) * (PE_nmas[i+tau_iter] - PE_nma_ave);
    }
    int num_iterations = iterations - tau_iter;
    double correlation = correlation_sum / PE_nma_var / num_iterations;
    
    sprintf(&correlation_buf[correlation_buf_offset],
	    "%+.5e     %+.5e\n", tau_iter*dt, correlation);
    correlation_buf_offset += 30;
    printf("correlation function progress: %d / %d, %g%%\r",
	   tau_iter, max_tau_iter, ((double) tau_iter) / max_tau_iter * 100);
    fflush(NULL);
  }
  
  FILE* correlation_data_file = fopen("pe_nma_correlation_function.txt", "w");
  fprintf(correlation_data_file, correlation_buf);
  fclose(correlation_data_file);
}
