#include "default_parameters.h"
#include "../dynein_struct.h"

void log_ts(void* dyn, State s, void* data, int iteration) {
  printf("another ts went by, it is now %d.\n", iteration);
  ((int*) data)[iteration] = iteration;
}

int main() {
  double runtime = dt*25;
  double test_position[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  void* data = malloc(ceil(runtime/dt) * sizeof(double));
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, test_position, log_ts, data);
  printf("data[5] = %d\n", ((int*) data)[5]);
  return EXIT_SUCCESS;
}
