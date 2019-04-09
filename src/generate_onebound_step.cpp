#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "default_parameters.h"
#include "dynein_struct.h"
#include "simulation_defaults.h"





int main(){
  double init_position[] = {120.0*M_PI / 180.0, 120.0*M_PI / 180.0, 0, 0, 8}; // nma_init, fma_init, nbx, nby, L

  simulateOneboundOnly(1.0, init_position); 
  return 0;
}



