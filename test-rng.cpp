#include "MersenneTwister.h"
#include <stdio.h>

int main() {
  double RAND_INIT_SEED = 0;
  MTRand* rand = new MTRand(RAND_INIT_SEED);
  long long count = 0;
  long long i = 0;
  long long iters = 1e14;
  while (i<iters) {
    i++;
    if (rand->rand() == 0) count++;
    if (i % (long long) 1e10 == 0) printf("count: %lld, iters: %lld, fraction: %g\n", count, i, double(count)/i);
  }
}
