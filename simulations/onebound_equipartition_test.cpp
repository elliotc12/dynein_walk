#include <cassert>
#include <limits>
#include <stdio.h>
#include <stdlib.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

void store_onebound_PEs(void* dyn, State s, void* job_msg, void* job_data, int iteration) {
  assert(s == NEARBOUND);
  Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
  ((double*) job_data)[4*iteration + 0] = dyn_ob->PE_bba;
  ((double*) job_data)[4*iteration + 1] = dyn_ob->PE_bma;
  ((double*) job_data)[4*iteration + 2] = dyn_ob->PE_ta;
  ((double*) job_data)[4*iteration + 3] = dyn_ob->PE_uma;
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
  void* data = malloc(iterations * sizeof(double) * 4);
  
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, test_position, store_onebound_PEs, NULL, data);

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

  printf("bba ratio: %g\n", get_average(bba_PEs, iterations)/(0.5*kb*T));
  printf("bma ratio: %g\n", get_average(bma_PEs, iterations)/(0.5*kb*T));
  printf("ta  ratio: %g\n", get_average(ta_PEs, iterations)/(0.5*kb*T));
  printf("uma ratio: %g\n", get_average(uma_PEs, iterations)/(0.5*kb*T));
}
