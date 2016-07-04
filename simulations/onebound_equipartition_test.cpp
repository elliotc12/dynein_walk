#include <cassert>
#include <limits>

#include "../default_parameters.h"
#include "../dynein_struct.h"
#include "simulations.h"

void get_average_PE(int iterations, double* averages) {
  double runtime = dt*iterations;
  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double test_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};
  void* data = malloc(iterations * sizeof(double) * 4);

  double* bba_PE_seed_ave = (double*) malloc(iterations * sizeof(double));
  double* bma_PE_seed_ave = (double*) malloc(iterations * sizeof(double));
  double*  ta_PE_seed_ave = (double*) malloc(iterations * sizeof(double));
  double* uma_PE_seed_ave = (double*) malloc(iterations * sizeof(double));
  
  for (int i = 0; i < iterations; i++) {
    bba_PE_seed_ave[i] = 0;
    bma_PE_seed_ave[i] = 0;
    ta_PE_seed_ave[i] = 0;
    uma_PE_seed_ave[i] = 0;
  }

  const int seeds[] = {0, 1, 2, 3, 4};
  int seed_len = sizeof(seeds) / sizeof(int);

  for (int r = 0; r < seed_len; r++) {
    RAND_INIT_SEED = seeds[r];
    simulate(runtime, RAND_INIT_SEED, NEARBOUND, test_position, store_onebound_PEs, &iterations, data);
  
    for (int i = 0; i < iterations; i++) {
      bba_PE_seed_ave[i] += ((double*) data)[4*i + 0];
      bma_PE_seed_ave[i] += ((double*) data)[4*i + 1];
      ta_PE_seed_ave[i]  += ((double*) data)[4*i + 2];
      uma_PE_seed_ave[i] += ((double*) data)[4*i + 3];
    }
  }

  for (int i = 0; i < iterations; i++) {
    bba_PE_seed_ave[i] /= seed_len;
    bma_PE_seed_ave[i] /= seed_len;
    ta_PE_seed_ave[i] /= seed_len;
    uma_PE_seed_ave[i] /= seed_len;
  }

  averages[0] = get_average(bba_PE_seed_ave, iterations);
  averages[1] = get_average(bma_PE_seed_ave, iterations);
  averages[2] = get_average(ta_PE_seed_ave, iterations);
  averages[3] = get_average(uma_PE_seed_ave, iterations);
}

int main() {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  double* low_T_PEs = (double*) malloc(4 * sizeof(double));
  double* high_T_PEs = (double*) malloc(4 * sizeof(double));

  int high_T = 10000;
  int low_T  = 100;
  
  int iterations = 1e7;
  T = low_T;
  get_average_PE(iterations, low_T_PEs);
  T = high_T;
  get_average_PE(iterations, high_T_PEs);

  double d_bba_PE = high_T_PEs[0] - low_T_PEs[0];
  double d_bma_PE = high_T_PEs[1] - low_T_PEs[1];
  double d_ta_PE =  high_T_PEs[2] - low_T_PEs[2];
  double d_uma_PE = high_T_PEs[3] - low_T_PEs[3];
  double d_T = high_T - low_T;

  printf("bba PE/0.5KBT ratio: %g\n", d_bba_PE / (0.5*kb*d_T));
  printf("bma PE/0.5KBT ratio: %g\n", d_bma_PE / (0.5*kb*d_T));
  printf("ta PE/0.5KBT ratio: %g\n", d_ta_PE / (0.5*kb*d_T));
  printf("uma PE/0.5KBT ratio: %g\n", d_uma_PE / (0.5*kb*d_T));
  printf("total PE/0.5KBT ratio: %g\n", (d_bba_PE + d_bma_PE + d_ta_PE + d_uma_PE) / (4*0.5*kb*d_T));
  
  return EXIT_SUCCESS;
}
