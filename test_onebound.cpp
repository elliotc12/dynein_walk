#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"

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

int main() {
  double bba_eq = bothbound_pre_powerstroke_internal_angles.bba;
  double ba_eq  = bothbound_pre_powerstroke_internal_angles.ba;
  double ta_eq  = bothbound_pre_powerstroke_internal_angles.ta;
  double fa_eq  = bothbound_pre_powerstroke_internal_angles.fa;

  forces no_forces =    {0,0,0,0,0,0,0,0,0,0}; // bbx, bby, bmx, bmy, ...
  forces right_forces = {1,0,1,0,1,0,1,0,1,0};
  forces left_forces =  {-1,0,-1,0,-1,0,-1,0,-1,0};
  forces up_forces =    {0,1,0,1,0,1,0,1,0,1};

  int num_failures = 0;
  {
    Dynein* dyn = new Dynein(bba_eq,                                 // starting bba
                             bba_eq + ba_eq - M_PI,                  // starting bma
                             bba_eq + ba_eq - M_PI + ta_eq,          // starting fma
                             bba_eq + ba_eq - M_PI + ta_eq + M_PI - fa_eq, //    fba
			     NEARBOUND,       // starting state
			     &no_forces,      // optional specified internal forces
			     &no_forces,      // optional specified brownian forces
		             (equilibrium_angles*) NULL);   // optional eq angles


    // Dynein in normal prepowerstroke conformation,
    // check if velocities agree with definitions.

    printf("Test: Dynein prepowerstroke conformation, "
	   "no internal forces, no Brownian forces.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
    num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
    num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
    num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
    num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);
    num_failures += test("Is d_fby zero", dyn->get_d_fby(), 0);

    num_failures += test("Is d_bba zero", dyn->get_d_bba(), 0);
    num_failures += test("Is d_bma zero", dyn->get_d_bma(), 0);
    num_failures += test("Is d_fma zero", dyn->get_d_fma(), 0);
    num_failures += test("Is d_fba zero", dyn->get_d_fba(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein((90.0 / 180) * M_PI,
                             (90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
			     NEARBOUND,
			     &no_forces,
			     &no_forces,
			     NULL);

    printf("\nTest: Dynein vertical conformation, no internal forces, no Brownian forces.\n");
    num_failures += test("Is fbx zero", dyn->get_fbx(), 0);
    num_failures += test("Is fmx zero", dyn->get_fmx(), 0);
    num_failures += test("Is tx zero", dyn->get_tx(), 0);
    num_failures += test("Is bmx zero", dyn->get_bmx(), 0);
    num_failures += test("Is bbx zero", dyn->get_bbx(), 0);
    
    num_failures += test("Is fby fully extended", dyn->get_fby(), 2*ls + 2*lt);
    num_failures += test("Is fmy fully extended", dyn->get_fmy(), 1*ls + 2*lt);
    num_failures += test("Is ty fully extended", dyn->get_ty(), 1*ls + 1*lt);
    num_failures += test("Is bmy fully extended", dyn->get_bmy(), 1*ls);
    num_failures += test("Is bby fully extended", dyn->get_bby(), 0);
    
    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
    num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
    num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
    num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
    num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);
    num_failures += test("Is d_fby zero", dyn->get_d_fby(), 0);

    num_failures += test("Is d_bba zero", dyn->get_d_bba(), 0);
    num_failures += test("Is d_bma zero", dyn->get_d_bma(), 0);
    num_failures += test("Is d_fma zero", dyn->get_d_fma(), 0);
    num_failures += test("Is d_fba zero", dyn->get_d_fba(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein((0.0 / 180) * M_PI,
                             (0.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
			     NEARBOUND,
			     &no_forces,
			     &no_forces,
			     NULL);
    
    printf("\nTest: Dynein horizontal conformation, no internal forces, no Brownian forces.\n");
    num_failures += test("Is fby zero", dyn->get_fby(), 0);
    num_failures += test("Is fmy zero", dyn->get_fmy(), 0);
    num_failures += test("Is ty zero", dyn->get_ty(), 0);
    num_failures += test("Is bmy zero", dyn->get_bmy(), 0);
    num_failures += test("Is bby zero", dyn->get_bby(), 0);

    num_failures += test("Is fbx fully extended", dyn->get_fbx(), 2*ls + 2*lt);
    num_failures += test("Is fmx fully extended", dyn->get_fmx(), 1*ls + 2*lt);
    num_failures += test("Is tx fully extended", dyn->get_tx(), 1*ls + 1*lt);
    num_failures += test("Is bmx fully extended", dyn->get_bmx(), 1*ls);
    num_failures += test("Is bbx fully extended", dyn->get_bbx(), 0);

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
    num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
    num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
    num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
    num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);

    num_failures += test("Is d_bba zero", dyn->get_d_bba(), 0);
    num_failures += test("Is d_bma zero", dyn->get_d_bma(), 0);
    num_failures += test("Is d_fma zero", dyn->get_d_fma(), 0);
    num_failures += test("Is d_fba zero", dyn->get_d_fba(), 0);
    
    free(dyn);
  }
  
  {
    Dynein* dyn = new Dynein(0.5*M_PI, 0.5*M_PI, 0.5*M_PI, 0.5*M_PI,
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL);

    printf("\nTest: Dynein vertical, no internal forces, Brownian forces in positive x direction.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test_noteq("Is d_bmx nonzero", dyn->get_d_bmx(), 0);
    num_failures += test_noteq("Is d_tx nonzero", dyn->get_d_tx(), 0);
    num_failures += test_noteq("Is d_fmx nonzero", dyn->get_d_fmx(), 0);
    num_failures += test_noteq("Is d_fbx nonzero", dyn->get_d_fbx(), 0);

    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
    num_failures += test("Is d_fby zero", dyn->get_d_fby(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein(0.5*M_PI, 0.5*M_PI, 0.5*M_PI, 0.5*M_PI,
			     NEARBOUND,
			     &no_forces,
			     &up_forces,
			     NULL);

    printf("\nTest: Dynein vertical, no internal forces, Brownian forces in positive y direction.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
    num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);

    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
    num_failures += test("Is d_fby zero", dyn->get_d_fby(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein(0, 0, 0, 0,
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL);

    printf("\nTest: Dynein horizontal, no internal forces, Brownian forces in positive x direction.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
    num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);

    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
    num_failures += test("Is d_fby zero", dyn->get_d_fby(), 0);

    free(dyn);
  }


  {
    Dynein* dyn = new Dynein(0, 0, 0, 0,
			     NEARBOUND,
			     &no_forces,
			     &up_forces,
			     NULL);

    printf("\nTest: Dynein horizontal, no internal forces, Brownian forces in positive y direction.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
    num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);

    num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
    num_failures += test_noteq("Is d_bmy nonzero", dyn->get_d_bmy(), 0);
    num_failures += test_noteq("Is d_ty nonzero", dyn->get_d_ty(), 0);
    num_failures += test_noteq("Is d_fmy nonzero", dyn->get_d_fmy(), 0);
    num_failures += test_noteq("Is d_fby nonzero", dyn->get_d_fby(), 0);

    free(dyn);
  }

  // { Test fails, cannot be in both leftbound and prepowerstroke states simultaneously
  //   Dynein* dyn = new Dynein(bba_eq,
  //                            bba_eq + ba_eq - M_PI,
  //                            bba_eq + ba_eq - M_PI + ta_eq,
  //                            bba_eq + ba_eq - M_PI + ta_eq + M_PI - fa_eq,
  //                            NEARBOUND,
  // 			     NULL,
  // 			     &no_forces,
  // 			     NULL);

  //   printf("\nTest: Dynein equilibrium, checking internal forces, no Brownian forces.\n");

  //   num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
  //   num_failures += test("Is d_bmx zero", dyn->get_d_bmx(), 0);
  //   num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
  //   num_failures += test("Is d_fmx zero", dyn->get_d_fmx(), 0);
  //   num_failures += test("Is d_fbx zero", dyn->get_d_fbx(), 0);

  //   num_failures += test("Is d_bby zero", dyn->get_d_bby(), 0);
  //   num_failures += test("Is d_bmy zero", dyn->get_d_bmy(), 0);
  //   num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
  //   num_failures += test("Is d_fmy zero", dyn->get_d_fmy(), 0);
  //   num_failures += test("Is d_fby zero", dyn->get_d_fby(), 0);

  //   num_failures += test("Is f_bbx zero", dyn->get_internal().bbx, 0);
  //   num_failures += test("Is f_bmx zero", dyn->get_internal().bmx, 0);
  //   num_failures += test("Is f_tx zero",  dyn->get_internal().tx, 0);
  //   num_failures += test("Is f_fmx zero", dyn->get_internal().fmx, 0);
  //   num_failures += test("Is f_fbx zero", dyn->get_internal().fbx, 0);

  //   num_failures += test("Is f_bby zero", dyn->get_internal().bby, 0);
  //   num_failures += test("Is f_bmy zero", dyn->get_internal().bmy, 0);
  //   num_failures += test("Is f_ty zero",  dyn->get_internal().ty, 0);
  //   num_failures += test("Is f_fmy zero", dyn->get_internal().fmy, 0);
  //   num_failures += test("Is f_fby zero", dyn->get_internal().fby, 0);

  //   free(dyn);
  // }

  {
    Dynein* dyn = new Dynein(bba_eq,
                             bba_eq + ba_eq - M_PI,
                             bba_eq + ba_eq - M_PI + ta_eq,
                             bba_eq + ba_eq - M_PI + ta_eq + M_PI - fa_eq,
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL);
    
    printf("\nTest: Dynein prepowerstroke conformation, no internal forces, Brownian forces in positive x direction.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test_noteq("Is d_bmx nonzero", dyn->get_d_bmx(), 0);
    num_failures += test_noteq("Is d_fmx nonzero", dyn->get_d_fmx(), 0);
    num_failures += test_noteq("Is d_fbx nonzero", dyn->get_d_fbx(), 0);
    
    free(dyn);
  }
  
  {
    Dynein* dyn = new Dynein(bba_eq,
                             bba_eq + ba_eq - M_PI,
                             bba_eq + ba_eq - M_PI + ta_eq,
                             bba_eq + ba_eq - M_PI + ta_eq + M_PI - fa_eq,
			     NEARBOUND,
			     &left_forces,
			     &no_forces,
			     NULL);
    
    printf("\nTest: Dynein prepowerstroke conformation, internal forces in negative x direction, no Brownian forces.\n");

    num_failures += test("Is d_bbx zero", dyn->get_d_bbx(), 0);
    num_failures += test_noteq("Is d_bmx nonzero", dyn->get_d_bmx(), 0);
    num_failures += test_noteq("Is d_fmx nonzero", dyn->get_d_fmx(), 0);
    num_failures += test_noteq("Is d_fbx nonzero", dyn->get_d_fbx(), 0);
    
    free(dyn);
  }
  
  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  exit(num_failures);
}
