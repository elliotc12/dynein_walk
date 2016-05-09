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
    return 1;
  } else {
    printf("%30s: FAIL! %g != %g.\n", msg, one, two);
    return 0;
  }
}

int test_noteq(const char *msg, float one, float two) {
  num_tests += 1;
  if (!equal(one, two)) {
    printf("%30s: pass, %g != %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%30s: FAIL! %g == %g.\n", msg, one, two);
    return 0;
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

    printf("Testing upwards line conformation\n");

    bothbound_equilibrium_angles line_eq_angles = {
      M_PI/2, M_PI, 0, M_PI, M_PI/2
    };

    Dynein_bothbound dyn_bb(M_PI,              // nma_init
                            M_PI,              // fma_init
                            0,                 // nbx_init
                            0,                 // nby_init
                            1e-25,             // L
                            NULL,              // internal forces
			    NULL,              // brownian forces
			    &line_eq_angles,   // equilibrium angles
			    rand);             // MTRand

    if (test("f.nbx zero?", dyn_bb.get_internal().nbx, 0) == 0) num_failures++;
    if (test("f.nmx zero?", dyn_bb.get_internal().nmx, 0) == 0) num_failures++;
    if (test("f.tx  zero?", dyn_bb.get_internal().tx , 0) == 0) num_failures++;
    if (test("f.fmx zero?", dyn_bb.get_internal().fmx, 0) == 0) num_failures++;
    if (test("f.fbx zero?", dyn_bb.get_internal().fbx, 0) == 0) num_failures++;
    if (test("f.nby zero?", dyn_bb.get_internal().nby, 0) == 0) num_failures++;
    if (test("f.nmy zero?", dyn_bb.get_internal().nmy, 0) == 0) num_failures++;
    if (test("f.ty  zero?", dyn_bb.get_internal().ty , 0) == 0) num_failures++;
    if (test("f.fmy zero?", dyn_bb.get_internal().fmy, 0) == 0) num_failures++;
    if (test("f.fby zero?", dyn_bb.get_internal().fby, 0) == 0) num_failures++;

    if (test("nbx zero?", dyn_bb.get_nbx(), 0) == 0) num_failures++;
    if (test("nmx zero?", dyn_bb.get_nmx(), 0) == 0) num_failures++;
    if (test("tx  zero?", dyn_bb.get_tx(), 0) == 0) num_failures++;
    if (test("fmx zero?", dyn_bb.get_fmx(), 0) == 0) num_failures++;
    if (test("fbx zero?", dyn_bb.get_fbx(), 0) == 0) num_failures++;

    if (test("nby zero?", dyn_bb.get_nby(), 0) == 0) num_failures++;
    if (test_noteq("nmy nonzero?", dyn_bb.get_nmy(), 0) == 0) num_failures++;
    if (test_noteq("ty  nonzero?", dyn_bb.get_ty(), 0) == 0) num_failures++;
    if (test_noteq("fmy nonzero?", dyn_bb.get_fmy(), 0) == 0) num_failures++;
    if (test("fby zero?", dyn_bb.get_fby(), 0) == 0) num_failures++;
  }

  { /*** Upwards line conformation ***/

    printf("Testing upwards line conformation\n");

    bothbound_equilibrium_angles natural_eq_angles = {
      108*M_PI/180,
      108*M_PI/180,
      108*M_PI/180,
      2*M_PI - 108*M_PI/180,
      M_PI - 108*M_PI/180
    };

    Dynein_bothbound dyn_bb(M_PI,              // nma_init
                            M_PI,              // fma_init
                            0,                 // nbx_init
                            0,                 // nby_init
                            1e-25,             // L
                            NULL,              // internal forces
			    NULL,              // brownian forces
			    &line_eq_angles,   // equilibrium angles
			    rand);             // MTRand

    if (test("f.nbx zero?", dyn_bb.get_internal().nbx, 0) == 0) num_failures++;
    if (test("f.nmx zero?", dyn_bb.get_internal().nmx, 0) == 0) num_failures++;
    if (test("f.tx  zero?", dyn_bb.get_internal().tx , 0) == 0) num_failures++;
    if (test("f.fmx zero?", dyn_bb.get_internal().fmx, 0) == 0) num_failures++;
    if (test("f.fbx zero?", dyn_bb.get_internal().fbx, 0) == 0) num_failures++;
    if (test("f.nby zero?", dyn_bb.get_internal().nby, 0) == 0) num_failures++;
    if (test("f.nmy zero?", dyn_bb.get_internal().nmy, 0) == 0) num_failures++;
    if (test("f.ty  zero?", dyn_bb.get_internal().ty , 0) == 0) num_failures++;
    if (test("f.fmy zero?", dyn_bb.get_internal().fmy, 0) == 0) num_failures++;
    if (test("f.fby zero?", dyn_bb.get_internal().fby, 0) == 0) num_failures++;

    if (test("nbx zero?", dyn_bb.get_nbx(), 0) == 0) num_failures++;
    if (test("nmx zero?", dyn_bb.get_nmx(), 0) == 0) num_failures++;
    if (test("tx  zero?", dyn_bb.get_tx(), 0) == 0) num_failures++;
    if (test("fmx zero?", dyn_bb.get_fmx(), 0) == 0) num_failures++;
    if (test("fbx zero?", dyn_bb.get_fbx(), 0) == 0) num_failures++;

    if (test("nby zero?", dyn_bb.get_nby(), 0) == 0) num_failures++;
    if (test_noteq("nmy nonzero?", dyn_bb.get_nmy(), 0) == 0) num_failures++;
    if (test_noteq("ty  nonzero?", dyn_bb.get_ty(), 0) == 0) num_failures++;
    if (test_noteq("fmy nonzero?", dyn_bb.get_fmy(), 0) == 0) num_failures++;
    if (test("fby zero?", dyn_bb.get_fby(), 0) == 0) num_failures++;
  }


  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  exit(num_failures);
}
