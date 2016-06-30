#include <stdio.h>

#include "../dynein_struct.h"
#include "../default_parameters.h"

int main main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  T = 100;
  int iterations = 1e8;

  int d_iter = 1e4;
  int max_tau_iter = floor(iterations * 0.2);
  
  //write_onebound_PE_correlation_function(iterations, d_iter, max_tau_iter);
  write_onebound_equipartition_ratio_tau_dependent(iterations, d_iter, max_tau_iter);
  return 0;
}
