#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "dynein_struct.h"

/*
 * Test dynein_bothobund state for errors in positions and force calculations.
 */

double EPSILON = 1e-12;
double runtime = 0; // dynein_struct.h needs this to be defined here

extern const double dt;

int equal(double v1, double v2) {
  if (fabs(v1) < EPSILON) { return fabs(v2) < EPSILON; }
  else return fabs(v1 - v2)/fabs(v1) <= EPSILON;
}

static int num_tests = 0;

int test(const char *msg, float one, float two) {
  num_tests += 1;
  if (equal(one, two)) {
    printf("%30s: pass, %g == %g.\n", msg, one, two);
    return 0;
  } else {
    printf("%30s: FAIL! %g != %g.\n", msg, one, two);
    return 1;
  }
}

int test_noteq(const char *msg, float one, float two) {
  num_tests += 1;
  if (!equal(one, two)) {
    printf("%30s: pass, %g != %g.\n", msg, one, two);
    return 0;
  } else {
    printf("%30s: FAIL! %g == %g.\n", msg, one, two);
    return 1;
  }
}

/* ******************************** MAIN **************************************** */

int main(int argvc, char **argv) {

  printf("****************Starting Bothbound Test****************\n\n");

  MTRand* rand = new MTRand(RAND_INIT_SEED);

  // double nba_eq  = bothbound_post_powerstroke_internal_angles.nba;
  // double nma_eq  = bothbound_post_powerstroke_internal_angles.nma;
  // double ta_eq   = bothbound_post_powerstroke_internal_angles.ta;
  // double fma_eq  = bothbound_post_powerstroke_internal_angles.fma;
  // double fba_eq  = bothbound_post_powerstroke_internal_angles.fba;

  int num_failures = 0;

  { /*** Upwards line conformation ***/

    bothbound_equilibrium_angles line_eq_angles = {
      M_PI/2, M_PI, 0, M_PI, M_PI/2
    };

    Dynein_bothbound dyn_bb(M_PI,              // nma_init
                            M_PI,              // fma_init
                            0,                 // nbx_init
                            0,                 // nby_init
                            0,                 // L
                            NULL,              // internal forces
			    NULL,              // brownian forces
			    &line_eq_angles,   // equilibrium angles
			    rand);             // MTRand

    if (test("Testing internal force.\n", dyn_bb.get_internal().nmx, 0))
      num_failures++;
  }

  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  exit(num_failures);
}
