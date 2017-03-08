#include "MersenneTwister.h"
#include <stdio.h>

int main() {
  double RAND_INIT_SEED = 0;
  MTRand* rand = new MTRand(RAND_INIT_SEED);
  double iters = 1e10;
  double count = 0;
  for (int i=0; i<iters; i++) {
    if (rand->rand() == 0) count++;
    if ((i % ((int) 1e8)) == 0) printf("count: %g, iters: %d, fraction: %g\n", count, i, count/iters);
  }
}
