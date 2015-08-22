#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

#define EPSILON 1e-12

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

int main() {
	forces zero_f;
	zero_f.r_blx = 0;     zero_f.r_bly = 0;
	zero_f.r_mlx = 0;     zero_f.r_mly = 0;
	zero_f.r_tx  = 0;     zero_f.r_ty  = 0;
	zero_f.r_mrx = 0;     zero_f.r_mry = 0;
	zero_f.r_brx = 0;     zero_f.r_bry = 0;

	zero_f.f_blx = 0;     zero_f.f_bly = 0;
	zero_f.f_mlx = 0;     zero_f.f_mly = 0;
	zero_f.f_tx  = 0;     zero_f.f_ty  = 0;
	zero_f.f_mrx = 0;     zero_f.f_mry = 0;
	zero_f.f_brx = 0;     zero_f.f_bry = 0;

	printf("\n\nRunning dynein_struct member value check...\n");

	// Dynein in normal pentagon conformation, do the velocities agree with the velocity definitions?
  int num_failures = 0;
  {
    Dynein* dyn = new Dynein((108.0 / 180) * M_PI,
                             (36.0 / 180) * M_PI,
                             -(144.0 / 180) * M_PI,
                             -(72.0 / 180) * M_PI,
                             zero_f);

    num_failures += test("Is d_blx defined properly", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly defined properly", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx defined properly", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly defined properly", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx defined properly", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry defined properly", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx defined properly", dyn->get_d_brx(), 0);
    num_failures += test("Is d_bry defined properly", dyn->get_d_bry(), 0);
  }

  {
    forces f = zero_f;
    f.f_bly = 1;
    f.f_mly = 4;
    f.f_ty  = 3;
    f.f_mry = 7;
    f.f_bry = 9.3;
    Dynein* dyn = new Dynein((90.0 / 180) * M_PI,
                             (90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
                             f);
    printf("\n");
    num_failures += test("it is standing up (y)", dyn->get_bry(), 4*10);
    num_failures += test("it is standing up (x)", dyn->get_brx(), 0);

    num_failures += test("Is d_blx defined properly", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly defined properly", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx defined properly", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly defined properly", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx defined properly", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry defined properly", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx defined properly", dyn->get_d_brx(), 0);
    num_failures += test("Is d_bry defined properly", dyn->get_d_bry(), 0);
  }

  {
    forces f = zero_f;
    f.f_bly = 0;
    f.f_mly = 0;
    f.f_ty  = 0;
    f.f_mry = 0;
    f.f_bry = 9.3;
    Dynein* dyn = new Dynein((0.0 / 180) * M_PI,
                             (0.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
                             f);
    printf("\n");
    num_failures += test("it is lying down (bry)", dyn->get_bry(), 0);
    num_failures += test("it is lying down (mry)", dyn->get_mry(), 0);
    num_failures += test("it is lying down (ty)", dyn->get_ty(), 0);
    num_failures += test("it is lying down (mly)", dyn->get_mly(), 0);
    num_failures += test("it is lying down (bly)", dyn->get_bly(), 0);

    num_failures += test("it is lying down (x)", dyn->get_brx(), 4*10);
    num_failures += test("it is lying down (x)", dyn->get_mrx(), 3*10);
    num_failures += test("it is lying down (x)", dyn->get_tx(), 2*10);
    num_failures += test("it is lying down (x)", dyn->get_mlx(), 1*10);
    num_failures += test("it is lying down (x)", dyn->get_blx(), 0);

    num_failures += test("Is d_blx defined properly", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly defined properly", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx defined properly", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly defined properly", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx defined properly", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry defined properly", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx defined properly", dyn->get_d_brx(), 0);
    num_failures += test_noteq("Is d_bry defined properly", dyn->get_d_bry(), 0);
    num_failures += test_noteq("Is d_bra defined properly", dyn->get_d_bra(), 0);
  }

	// TODO: test other weird states
	// Vertical dynein with only f_brx nonzero: all velocities zero except bra
	// Horizontal dynein with only f_bry nonzero: all velocities zero except bra

  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
	exit(num_failures);
}













