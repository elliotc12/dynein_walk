#include <cassert>
#include <limits>

#include "../default_parameters.h"
#include "../dynein_struct.h"
#include "simulations.h"

int main() {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();

  double* low_T_PEs = (double*) malloc(4 * sizeof(double));
  double* high_T_PEs = (double*) malloc(4 * sizeof(double));

  int high_T = 10000;
  int low_T  = 100;
  
  int iterations = 1e7;
  const int seeds[] = {0, 1, 2, 3, 4};
  int seed_len = sizeof(seeds) / sizeof(int);
  
  T = low_T;
  get_average_PE(iterations, low_T_PEs, seeds, seed_len);
  
  T = high_T;
  get_average_PE(iterations, high_T_PEs, seeds, seed_len);

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
