#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"
#include "simulations.h"

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;
  int iterations = 1e6;

  int d_iter = 1e4;
  int max_tau_iter = floor(iterations * 0.5);
  int min_tau_iter = floor(iterations * 0.1);
  
  write_onebound_PE_correlation_function(iterations, d_iter, max_tau_iter);
  write_onebound_equipartition_ratio_per_tau(10*iterations, 10*d_iter, 10*min_tau_iter, 10*max_tau_iter);
  return 0;
}
