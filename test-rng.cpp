#include "MersenneTwister.h"
#include <stdio.h>

int main() {
  double RAND_INIT_SEED = 0;
  MTRand* rand = new MTRand(RAND_INIT_SEED);
  FILE* fd = fopen("rnglog.txt", "w");
  long long count = 0;
  long long i = 0;
  while (true) {
    i++;
    if (rand->rand() == 0) count++;
    if ((i % ((int) 1e8)) == 0) fprintf(fd, "count: %lld, iters: %lld, fraction: %g\n", count, i, count/i);
  }
}
