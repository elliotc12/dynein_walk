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

int equal(double v1, double v2, double epsilon = EPSILON) {
  if (fabs(v1) < epsilon) { return fabs(v2) < epsilon; }
  else return fabs(v1 - v2)/fabs(v1) <= epsilon;
}

static int num_tests = 0;

int test(const char *msg, float one, float two, double epsilon = EPSILON) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%30s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (equal(one, two, epsilon)) {
    printf("%30s: pass, %g == %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%30s: FAIL! %g != %g.\n", msg, one, two);
    return 0;
  }
}

int test_noteq(const char *msg, float one, float two) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%30s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (!equal(one, two)) {
    printf("%30s: pass, %g != %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%30s: FAIL! %g == %g.\n", msg, one, two);
    return 0;
  }
}

int test_greater(const char *msg, float one, float two, double epsilon = EPSILON) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%30s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (one > two) {
    printf("%30s: pass, %g > %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%30s: FAIL! %g <= %g.\n", msg, one, two);
    return 0;
  }
}

int test_less(const char *msg, float one, float two, double epsilon = EPSILON) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%30s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (one < two) {
    printf("%30s: pass, %g < %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%30s: FAIL! %g >= %g.\n", msg, one, two);
    return 0;
  }
}

/* ******************************** MAIN **************************************** */

int main(int argvc, char **argv) {

  const double pN = 1e-12; // a pico-Newton is a reasonable amount of force

  printf("****************Starting Bothbound Test****************\n\n");

  MTRand* rand = new MTRand(RAND_INIT_SEED);

  int num_failures = 0;

  bothbound_forces no_forces = {0,0,0,0,0,0,0,0,0,0}; // to eliminate brownian/equilibrium forces

  { printf("**Upwards line conformation with no forces**\n");

    bothbound_equilibrium_angles line_eq_angles = {
      M_PI/2, M_PI, 0, M_PI, M_PI/2
    };

    Dynein_bothbound dyn_bb(M_PI - 1e-15,      // nma_init
                            M_PI - 1e-15,      // fma_init
                            0,                 // nbx_init
                            0,                 // nby_init
                            1e-25,             // L
                            NULL,              // internal forces
			    &no_forces,        // brownian forces
			    &line_eq_angles,   // equilibrium angles
			    rand);             // MTRand

    printf("\tTesting coordinate system:\n");
    if (!test("nbx zero?", dyn_bb.get_nbx(), 0)) num_failures++;
    if (!test("nmx zero?", dyn_bb.get_nmx(), 0)) num_failures++;
    if (!test("tx  zero?", dyn_bb.get_tx(), 0)) num_failures++;
    if (!test("fmx zero?", dyn_bb.get_fmx(), 0)) num_failures++;
    if (!test("fbx zero?", dyn_bb.get_fbx(), 0)) num_failures++;

    if (!test("nby zero?", dyn_bb.get_nby(), 0)) num_failures++;
    if (!test_noteq("nmy nonzero?", dyn_bb.get_nmy(), 0)) num_failures++;
    if (!test_noteq("ty  nonzero?", dyn_bb.get_ty(), 0)) num_failures++;
    if (!test_noteq("fmy nonzero?", dyn_bb.get_fmy(), 0)) num_failures++;
    if (!test("fby zero?", dyn_bb.get_fby(), 0)) num_failures++;

    printf("\tTesting force definitions:\n");
    if (!test("f.nbx zero?", dyn_bb.get_internal().nbx, 0)) num_failures++;
    if (!test("f.nmx zero?", dyn_bb.get_internal().nmx, 0)) num_failures++;
    if (!test("f.tx  zero?", dyn_bb.get_internal().tx , 0)) num_failures++;
    if (!test("f.fmx zero?", dyn_bb.get_internal().fmx, 0)) num_failures++;
    if (!test("f.fbx zero?", dyn_bb.get_internal().fbx, 0)) num_failures++;
    if (!test("f.nby zero?", dyn_bb.get_internal().nby, 0)) num_failures++;
    if (!test("f.nmy zero?", dyn_bb.get_internal().nmy, 0)) num_failures++;
    if (!test("f.ty  zero?", dyn_bb.get_internal().ty , 0)) num_failures++;
    if (!test("f.fmy zero?", dyn_bb.get_internal().fmy, 0)) num_failures++;
    if (!test("f.fby zero?", dyn_bb.get_internal().fby, 0)) num_failures++;
  }

   { printf("\n**House conformation with equilateral roof, outwards forces**\n");

     bothbound_forces out_forces = {-pN*1000,0,-pN*1000,0,0,0,pN*1000,0,pN*1000,0};

     Dynein_bothbound dyn_bb(5*M_PI/6,     // nma_init
                            7*M_PI/6,      // fma_init
                            0,             // nbx_init
                            0,             // nby_init
                            Lt,            // L -- equilateral roof
                            &no_forces,    // internal forces
  			    &out_forces,   // brownian forces
  			    NULL,          // equilibrium angles
  			    rand);         // MTRand

     if (!test("nmx zero?", dyn_bb.get_nmx(), 0)) num_failures++;
     if (!test("fmx Lt?", dyn_bb.get_fmx(), Ls)) num_failures++;
     if (!test("nmy Ls?", dyn_bb.get_nmy(), Ls)) num_failures++;
     if (!test("fmy Ls?", dyn_bb.get_fmy(), Ls)) num_failures++;
     if (!test("tx Lt/2?", dyn_bb.get_tx(), Lt/2)) num_failures++;
     if (!test("ty Ls + sqrt(3)/2*Lt?", dyn_bb.get_ty(), Ls + sqrt(3)/2*Lt)) num_failures++;

     if (!test_less("d_nmx < 0?", dyn_bb.get_d_nmx(), 0)) num_failures++;
     if (!test("d_tx zero?", dyn_bb.get_d_tx(), 0)) num_failures++;
     if (!test_less("d_ty < 0?", dyn_bb.get_d_ty(), 0)) num_failures++;
     if (!test_greater("d_fmx > 0?", dyn_bb.get_d_fmx(), 0)) num_failures++;
  }

  { printf("\n**'Almost' upwards line conformation with +x forces**\n");

    bothbound_equilibrium_angles line_eq_angles = {
      M_PI/2, M_PI, 0, M_PI, M_PI/2
    };

    bothbound_forces x_forces =    {pN,0,pN,0,0,0,pN,0,pN,0}; // bbx, bby, bmx, bmy, ...

    Dynein_bothbound dyn_bb(M_PI - 1e-7,      // nma_init
                            M_PI - 1e-7,      // fma_init
                            0,                // nbx_init
                            0,                // nby_init
                            1e-25,            // L
                            NULL,             // internal forces
			    &x_forces,        // brownian forces
			    &line_eq_angles,  // equilibrium angles
			    rand);            // MTRand

    printf("\tTesting motor velocities:\n");
    if (!test_noteq("d_nmx_dt nonzero?", dyn_bb.get_d_nmx(), 0)) num_failures++;
    if (!test_noteq("d_tx_dt  nonzero?", dyn_bb.get_d_tx(), 0)) num_failures++;
    if (!test_noteq("d_fmx_dt nonzero?", dyn_bb.get_d_fmx(), 0)) num_failures++;

    if (!test("d_nmy_dt almost zero?", dyn_bb.get_d_nmy(), 0, 1e-4)) num_failures++;
    if (!test("d_ty_dt  almost zero?", dyn_bb.get_d_ty(), 0, 1e-4)) num_failures++;
    if (!test("d_fmy_dt almost zero?", dyn_bb.get_d_fmy(), 0, 1e-4)) num_failures++;
  }

  { printf("\n**Two tables with near/far domains flipped**\n");

    bothbound_equilibrium_angles left_table_eq_angles = {
      M_PI/2, M_PI/2, M_PI, 3*M_PI/2, M_PI/2
    };

    bothbound_equilibrium_angles right_table_eq_angles = {
      M_PI/2, 3*M_PI/2, -M_PI, M_PI/2, M_PI/2
    };

    Dynein_bothbound left_dyn_bb(M_PI/2,                 // nma_init
				 3*M_PI/2,               // fma_init
				 0,                      // nbx_init
				 0,                      // nby_init
				 2*Lt,                   // L
				 NULL,                   // internal forces
				 &no_forces,             // brownian forces
				 &left_table_eq_angles,  // equilibrium angles
				 rand);                  // MTRand

    Dynein_bothbound right_dyn_bb(3*M_PI/2,               // nma_init
				 M_PI/2,                  // fma_init
				 2*Lt,                    // nbx_init
				 0,                       // nby_init
				 -2*Lt,                   // L
				 NULL,                    // internal forces
				 &no_forces,              // brownian forces
				 &right_table_eq_angles,  // equilibrium angles
				 rand);                   // MTRand

    if (!test("left bx coords equal?",
	      left_dyn_bb.get_nbx(), right_dyn_bb.get_fbx())) num_failures++;
    if (!test("left mx coords equal?",
	      left_dyn_bb.get_nmx(), right_dyn_bb.get_fmx())) num_failures++;
    if (!test("tx coords equal?",
	      left_dyn_bb.get_tx(), right_dyn_bb.get_tx())) num_failures++;
    if (!test("right mx coords equal?",
	      left_dyn_bb.get_fmx(), right_dyn_bb.get_nmx())) num_failures++;
    if (!test("right bx coords equal?",
	      left_dyn_bb.get_fbx(), right_dyn_bb.get_nbx())) num_failures++;

    if (!test("left my coords equal?",
	      left_dyn_bb.get_nmy(), right_dyn_bb.get_fmy())) num_failures++;
    if (!test("ty coords equal?",
	      left_dyn_bb.get_ty(), right_dyn_bb.get_ty())) num_failures++;
    if (!test("right my coords equal?",
	      left_dyn_bb.get_fmy(), right_dyn_bb.get_nmy())) num_failures++;

    if (!test("left angle velocities equal?",
	      left_dyn_bb.get_d_nma(), right_dyn_bb.get_d_fma())) num_failures++;
    if (!test("right angle velocities equal?",
	      left_dyn_bb.get_d_fma(), right_dyn_bb.get_d_nma())) num_failures++;

    if (!test("left mx velocities equal?",
	      left_dyn_bb.get_d_nmx(), right_dyn_bb.get_d_fmx())) num_failures++;
    if (!test("tx velocities equal?",
	      left_dyn_bb.get_d_tx(), right_dyn_bb.get_d_tx())) num_failures++;
    if (!test("right mx velocities equal?",
	      left_dyn_bb.get_d_fmx(), right_dyn_bb.get_d_nmx())) num_failures++;

    if (!test("left my velocities equal?",
	      left_dyn_bb.get_d_nmy(), right_dyn_bb.get_d_fmy())) num_failures++;
    if (!test("ty velocities equal?",
	      left_dyn_bb.get_d_ty(), right_dyn_bb.get_d_ty())) num_failures++;
    if (!test("right my velocities equal?",
	      left_dyn_bb.get_d_fmy(), right_dyn_bb.get_d_nmy())) num_failures++;
  }

  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  exit(num_failures);
}
