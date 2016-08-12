#include <cassert>
#include <limits>

#include "../default_parameters.h"
#include "../dynein_struct.h"

void log_bothbound_PEs(void* dyn, State s, void* job_msg, void* job_data, int iteration) {
  assert(s == BOTHBOUND);
  Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
  ((double*) job_data)[5*iteration + 0] = dyn_bb->PE_nba;
  ((double*) job_data)[5*iteration + 1] = dyn_bb->PE_nma;
  ((double*) job_data)[5*iteration + 2] = dyn_bb->PE_ta;
  ((double*) job_data)[5*iteration + 3] = dyn_bb->PE_fma;
  ((double*) job_data)[5*iteration + 4] = dyn_bb->PE_fba;
}

int main() {
  T = 100;
  BOTHBOUND_UNBINDING_FORCE = std::numeric_limits<double>::infinity();
  
  int iterations = 1e7;
  double runtime = dt*iterations;
  bothbound_equilibrium_angles eq = bothbound_pre_powerstroke_internal_angles;
  double test_position[] = {eq.nma, eq.fma, 0.0, 0.0, Ls};
  void* data = malloc(iterations * sizeof(double) * 5);

  simulate(runtime, RAND_INIT_SEED, BOTHBOUND, test_position, log_bothbound_PEs, NULL, data);

  double* nba_PEs = (double*) malloc(iterations * sizeof(double));
  double* nma_PEs = (double*) malloc(iterations * sizeof(double));
  double*  ta_PEs = (double*) malloc(iterations * sizeof(double));
  double* fma_PEs = (double*) malloc(iterations * sizeof(double));
  double* fba_PEs = (double*) malloc(iterations * sizeof(double));
  
  for (int i = 0; i < iterations; i++) {
    nba_PEs[i] = ((double*) data)[5*i + 0];
    nma_PEs[i] = ((double*) data)[5*i + 1];
    ta_PEs[i]  = ((double*) data)[5*i + 2];
    fma_PEs[i] = ((double*) data)[5*i + 3];
    fba_PEs[i] = ((double*) data)[5*i + 4];
  }

  printf("0.*kb*T: %g\n", 0.5*kb*T);
  printf("nba_PE_ave: %g, ratio: %g\n", get_average(nba_PEs, iterations),
	 get_average(nba_PEs, iterations)/(0.5*kb*T));
  printf("nma_PE_ave: %g, ratio: %g\n", get_average(nma_PEs, iterations),
	 get_average(nma_PEs, iterations)/(0.5*kb*T));
  printf(" ta_PE_ave: %g, ratio: %g\n", get_average( ta_PEs, iterations),
	 get_average( ta_PEs, iterations)/(0.5*kb*T));
  printf("fma_PE_ave: %g, ratio: %g\n", get_average(fma_PEs, iterations),
	 get_average(fma_PEs, iterations)/(0.5*kb*T));
  printf("fba_PE_ave: %g, ratio: %g\n", get_average(fba_PEs, iterations),
	 get_average(fba_PEs, iterations)/(0.5*kb*T));
  
  return EXIT_SUCCESS;
}
