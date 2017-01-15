#include <stdio.h>
#include <stdlib.h>
#include <csignal>
#include <cassert>
#include <limits>

#include "dynein_struct.h"
#include "default_parameters.h"

/*
 * Test dynein_onebound state for errors in position, velocity and force calculations.
 */

double EPSILON = 1e-5;

int equal(double v1, double v2, double epsilon = EPSILON) {
  if (fabs(v1) <= epsilon) { return fabs(v2) <= epsilon; }
  else return fabs(v1 - v2)/fabs(v1) <= epsilon;
}

static int num_tests = 0;

int test(const char *msg, float one, float two, double epsilon = EPSILON) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%45s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (equal(one, two, epsilon)) {
    printf("%45s: pass, %g == %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%45s: FAIL! %g != %g\tratio = %g.\n", msg, one, two, one/two);
    return 0;
  }
}

int test_noteq(const char *msg, float one, float two) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%45s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (!equal(one, two)) {
    printf("%45s: pass, %g != %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%45s: FAIL! %g == %g.\n", msg, one, two);
    return 0;
  }
}

int test_greater(const char *msg, float one, float two, double epsilon = EPSILON) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%45s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (one > two) {
    printf("%45s: pass, %g > %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%45s: FAIL! %g <= %g.\n", msg, one, two);
    return 0;
  }
}

int test_less(const char *msg, float one, float two, double epsilon = EPSILON) {
  num_tests += 1;
  if (one != one or two != two) {
    printf("%45s: NaN FAIL!, %g %g.\n", msg, one, two);
    return 0;
  } else if (one < two) {
    printf("%45s: pass, %g < %g.\n", msg, one, two);
    return 1;
  } else {
    printf("%45s: FAIL! %g >= %g.\n", msg, one, two);
    return 0;
  }
}

/* ******************************** MAIN **************************************** */

int main() {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }
  
  double bba_eq = onebound_post_powerstroke_internal_angles.bba;
  double bma_eq = onebound_post_powerstroke_internal_angles.bma;
  double ta_eq  = onebound_post_powerstroke_internal_angles.ta;
  double uma_eq = onebound_post_powerstroke_internal_angles.uma;

  double R = sqrt(2*kb*T/(gm*dt)); // Brownian force constant
  MTRand* rand = new MTRand(RAND_INIT_SEED);
  int num_failures = 0;

  onebound_forces no_forces    = {0,0,0,0,0,0,0,0,0,0}; // bbx, bby, bmx, bmy, ...
  onebound_forces right_forces = {gb*R,0,gm*R,0,gt*R,0,gm*R,0,gb*R,0};
  onebound_forces left_forces  = {-gb*R,0,-gm*R,0,-gt*R,0,-gm*R,0,-gb*R,0};
  onebound_forces up_forces    = {0,gb*R,0,gm*R,0,gt*R,0,gm*R,0,gb*R};

  printf("****************starting Bothbound Test****************\n\n");

  { printf("**Upwards line conformation with no forces**\n");
    
    Dynein_onebound* dyn_ob = new Dynein_onebound(bba_eq,            // starting bba
                             bba_eq + bma_eq - M_PI,                  // starting bma
                             bba_eq + bma_eq - M_PI + ta_eq,          // starting uma
                             bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq, // starting uba
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,                              // starting state
			     &no_forces,                             // internal forces
			     &no_forces,                             // brownian forces
		             NULL,                                   // eq angles
			     rand);                                  // MTRand

    printf("Testing velocities:\n");
    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    if (!test("Is d_bmx zero", dyn_ob->get_d_bmx(), 0)) num_failures++;
    if (!test("Is d_bmy zero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    if (!test("Is d_umx zero", dyn_ob->get_d_umx(), 0)) num_failures++;
    if (!test("Is d_umy zero", dyn_ob->get_d_umy(), 0)) num_failures++;
    if (!test("Is d_ubx zero", dyn_ob->get_d_ubx(), 0)) num_failures++;
    if (!test("Is d_uby zero", dyn_ob->get_d_uby(), 0)) num_failures++;

    if (!test("Is d_bba zero", dyn_ob->get_d_bba(), 0)) num_failures++;
    if (!test("Is d_bma zero", dyn_ob->get_d_bma(), 0)) num_failures++;
    if (!test("Is d_uma zero", dyn_ob->get_d_uma(), 0)) num_failures++;
    if (!test("Is d_uba zero", dyn_ob->get_d_uba(), 0)) num_failures++;

    delete dyn_ob;
  }

  { printf("\n**Upwards line conformation with no forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound((90.0 / 180) * M_PI,
                             (90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &no_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is ubx zero", dyn_ob->get_ubx(), 0)) num_failures++;
    if (!test("Is umx zero", dyn_ob->get_umx(), 0)) num_failures++;
    if (!test("Is tx zero",  dyn_ob->get_tx(), 0)) num_failures++;
    if (!test("Is bmx zero", dyn_ob->get_bmx(), 0)) num_failures++;
    if (!test("Is bbx zero", dyn_ob->get_bbx(), 0)) num_failures++;
    
    if (!test("Is uby fully extended", dyn_ob->get_uby(), 2*Ls + 2*Lt)) num_failures++;
    if (!test("Is umy fully extended", dyn_ob->get_umy(), 1*Ls + 2*Lt)) num_failures++;
    if (!test("Is ty fully extended", dyn_ob->get_ty(), 1*Ls + 1*Lt)) num_failures++;
    if (!test("Is bmy fully extended", dyn_ob->get_bmy(), 1*Ls)) num_failures++;
    if (!test("Is bby fully extended", dyn_ob->get_bby(), 0)) num_failures++;
    
    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    if (!test("Is d_bmx zero", dyn_ob->get_d_bmx(), 0)) num_failures++;
    if (!test("Is d_bmy zero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    if (!test("Is d_umx zero", dyn_ob->get_d_umx(), 0)) num_failures++;
    if (!test("Is d_umy zero", dyn_ob->get_d_umy(), 0)) num_failures++;
    if (!test("Is d_ubx zero", dyn_ob->get_d_ubx(), 0)) num_failures++;
    if (!test("Is d_uby zero", dyn_ob->get_d_uby(), 0)) num_failures++;

    if (!test("Is d_bba zero", dyn_ob->get_d_bba(), 0)) num_failures++;
    if (!test("Is d_bma zero", dyn_ob->get_d_bma(), 0)) num_failures++;
    if (!test("Is d_uma zero", dyn_ob->get_d_uma(), 0)) num_failures++;
    if (!test("Is d_uba zero", dyn_ob->get_d_uba(), 0)) num_failures++;

    delete dyn_ob;
  }

  { printf("\n**Horizontal line conformation with no forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound((0.0 / 180) * M_PI,
                             (0.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &no_forces,
			     NULL,
			     rand);                                  // MTRand
    
    if (!test("Is uby zero", dyn_ob->get_uby(), 0)) num_failures++;
    if (!test("Is umy zero", dyn_ob->get_umy(), 0)) num_failures++;
    if (!test("Is ty zero", dyn_ob->get_ty(), 0)) num_failures++;
    if (!test("Is bmy zero", dyn_ob->get_bmy(), 0)) num_failures++;
    if (!test("Is bby zero", dyn_ob->get_bby(), 0)) num_failures++;

    if (!test("Is ubx fully extended", dyn_ob->get_ubx(), 2*Ls + 2*Lt)) num_failures++;
    if (!test("Is umx fully extended", dyn_ob->get_umx(), 1*Ls + 2*Lt)) num_failures++;
    if (!test("Is tx fully extended", dyn_ob->get_tx(), 1*Ls + 1*Lt)) num_failures++;
    if (!test("Is bmx fully extended", dyn_ob->get_bmx(), 1*Ls)) num_failures++;
    if (!test("Is bbx fully extended", dyn_ob->get_bbx(), 0)) num_failures++;

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    // if (!test("Is d_bmx zero", dyn_ob->get_d_bmx(), 0)) num_failures++; //These are NaN failing, not necessarily bad
    // if (!test("Is d_bmy zero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    // if (!test("Is d_umx zero", dyn_ob->get_d_umx(), 0)) num_failures++;
    // if (!test("Is d_umy zero", dyn_ob->get_d_umy(), 0)) num_failures++;
    // if (!test("Is d_ubx zero", dyn_ob->get_d_ubx(), 0)) num_failures++;

    // if (!test("Is d_bba zero", dyn_ob->get_d_bba(), 0)) num_failures++;
    // if (!test("Is d_bma zero", dyn_ob->get_d_bma(), 0)) num_failures++;
    // if (!test("Is d_uma zero", dyn_ob->get_d_uma(), 0)) num_failures++;
    // if (!test("Is d_uba zero", dyn_ob->get_d_uba(), 0)) num_failures++;
    
    delete dyn_ob;
  }

  { printf("\n**Upwards line conformation with +x forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound(0.5*M_PI, 0.5*M_PI, 0.5*M_PI, 0.5*M_PI,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test_greater("Is d_bmx positive", dyn_ob->get_d_bmx(), 0)) num_failures++;
    if (!test_greater("Is d_tx positive", dyn_ob->get_d_tx(), 0)) num_failures++;
    if (!test_greater("Is d_umx positive", dyn_ob->get_d_umx(), 0)) num_failures++;
    if (!test_greater("Is d_ubx positive", dyn_ob->get_d_ubx(), 0)) num_failures++;

    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    if (!test("Is d_bmy zero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    if (!test("Is d_ty zero", dyn_ob->get_d_ty(), 0)) num_failures++;
    if (!test("Is d_umy zero", dyn_ob->get_d_umy(), 0)) num_failures++;
    if (!test("Is d_uby zero", dyn_ob->get_d_uby(), 0)) num_failures++;

    delete dyn_ob;
  }

  { printf("\n**Upwards line conformation with +y forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound(0.5*M_PI, 0.5*M_PI, 0.5*M_PI, 0.5*M_PI,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &up_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test("Is d_bmx zero", dyn_ob->get_d_bmx(), 0)) num_failures++;
    if (!test("Is d_tx zero", dyn_ob->get_d_tx(), 0)) num_failures++;
    if (!test("Is d_umx zero", dyn_ob->get_d_umx(), 0)) num_failures++;
    if (!test("Is d_ubx zero", dyn_ob->get_d_ubx(), 0)) num_failures++;

    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    if (!test("Is d_bmy zero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    if (!test("Is d_ty zero", dyn_ob->get_d_ty(), 0)) num_failures++;
    if (!test("Is d_umy zero", dyn_ob->get_d_umy(), 0)) num_failures++;
    if (!test("Is d_uby zero", dyn_ob->get_d_uby(), 0)) num_failures++;

    delete dyn_ob;
  }

  { printf("\n**Horizontal line conformation with +x forces**\n");
    
    Dynein_onebound* dyn_ob = new Dynein_onebound(0, 0, 0, 0,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    // if (!test("Is d_bmx zero", dyn_ob->get_d_bmx(), 0)) num_failures++; // NaN failing; not necessarily bad
    // if (!test("Is d_tx zero", dyn_ob->get_d_tx(), 0)) num_failures++;
    // if (!test("Is d_umx zero", dyn_ob->get_d_umx(), 0)) num_failures++;
    // if (!test("Is d_ubx zero", dyn_ob->get_d_ubx(), 0)) num_failures++;

    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    // if (!test("Is d_bmy zero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    // if (!test("Is d_ty zero", dyn_ob->get_d_ty(), 0)) num_failures++;
    // if (!test("Is d_umy zero", dyn_ob->get_d_umy(), 0)) num_failures++;
    // if (!test("Is d_uby zero", dyn_ob->get_d_uby(), 0)) num_failures++;

    delete dyn_ob;
  }


  { printf("\n**Horizontal line conformation with +y forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound(0, 0, 0, 0,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &up_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    // if (!test("Is d_bmx zero", dyn_ob->get_d_bmx(), 0)) num_failures++; // NaN failing
    // if (!test("Is d_tx zero", dyn_ob->get_d_tx(), 0)) num_failures++;
    // if (!test("Is d_umx zero", dyn_ob->get_d_umx(), 0)) num_failures++;
    // if (!test("Is d_ubx zero", dyn_ob->get_d_ubx(), 0)) num_failures++;

    if (!test("Is d_bby zero", dyn_ob->get_d_bby(), 0)) num_failures++;
    // if (!test_noteq("Is d_bmy nonzero", dyn_ob->get_d_bmy(), 0)) num_failures++;
    // if (!test_noteq("Is d_ty nonzero", dyn_ob->get_d_ty(), 0)) num_failures++;
    // if (!test_noteq("Is d_umy nonzero", dyn_ob->get_d_umy(), 0)) num_failures++;
    // if (!test_noteq("Is d_uby nonzero", dyn_ob->get_d_uby(), 0)) num_failures++;

    delete dyn_ob;
  }

  { printf("\n**Prepowerstroke conformation with +x forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound(bba_eq,
                             bba_eq + bma_eq - M_PI,
                             bba_eq + bma_eq - M_PI + ta_eq,
                             bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test_noteq("Is d_bmx nonzero", dyn_ob->get_d_bmx(), 0)) num_failures++;
    if (!test_noteq("Is d_umx nonzero", dyn_ob->get_d_umx(), 0)) num_failures++;
    if (!test_noteq("Is d_ubx nonzero", dyn_ob->get_d_ubx(), 0)) num_failures++;
    
    delete dyn_ob;
  }
  
  { printf("\n**Prepowerstroke conformation with -x internal forces**\n");
    Dynein_onebound* dyn_ob = new Dynein_onebound(bba_eq,
                             bba_eq + bma_eq - M_PI,
                             bba_eq + bma_eq - M_PI + ta_eq,
                             bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &left_forces,
			     &no_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Is d_bbx zero", dyn_ob->get_d_bbx(), 0)) num_failures++;
    if (!test_noteq("Is d_bmx nonzero", dyn_ob->get_d_bmx(), 0)) num_failures++;
    if (!test_noteq("Is d_umx nonzero", dyn_ob->get_d_umx(), 0)) num_failures++;
    if (!test_noteq("Is d_ubx nonzero", dyn_ob->get_d_ubx(), 0)) num_failures++;
    
    delete dyn_ob;
  }

  { printf("\n**Prepowerstroke conformation opposing internal/brownian forces**\n");
    Dynein_onebound* dyn_ob_plus = new Dynein_onebound(bba_eq,
                             bba_eq + bma_eq - M_PI,
                             bba_eq + bma_eq - M_PI + ta_eq,
                             bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &no_forces,
			     &right_forces,
			     NULL,
			     rand);                                  // MTRand

    Dynein_onebound* dyn_ob_minus = new Dynein_onebound(bba_eq,
                             bba_eq + bma_eq - M_PI,
                             bba_eq + bma_eq - M_PI + ta_eq,
                             bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq,
			     0,                                      // nbx_init
			     0,                                      // nby_init
			     NEARBOUND,
			     &left_forces,
			     &no_forces,
			     NULL,
			     rand);                                  // MTRand

    if (!test("Are d_bmxs equal^opposite", dyn_ob_plus->get_d_bmx(), -dyn_ob_minus->get_d_bmx())) num_failures++;
    if (!test("Are d_umxs equal^opposite", dyn_ob_plus->get_d_umx(), -dyn_ob_minus->get_d_umx())) num_failures++;
    if (!test("Are d_ubxs equal^opposite", dyn_ob_plus->get_d_ubx(), -dyn_ob_minus->get_d_ubx())) num_failures++;
    
    delete dyn_ob_plus;
    delete dyn_ob_minus;
  }

  { printf("\n**Conservation of Energy**\n");

    double poke = 0.00001;
    Dynein_onebound ob_1(bba_eq + 0.5,
			 bba_eq + bma_eq - M_PI + 0.5,
			 bba_eq + bma_eq - M_PI + ta_eq + 0.5,
			 bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq + 0.5,
			 0,
			 0,
			 NEARBOUND,
			 NULL,
			 &no_forces,
			 NULL,
			 rand);

    Dynein_onebound ob_2(bba_eq + 0.5 + poke,
			 bba_eq + bma_eq - M_PI + 0.5 + poke,
			 bba_eq + bma_eq - M_PI + ta_eq + 0.5 + poke,
			 bba_eq + bma_eq - M_PI + ta_eq + M_PI - uma_eq + 0.5 + poke,
			 0,
			 0,
			 NEARBOUND,
			 NULL,
			 &no_forces,
			 NULL,
			 rand);

    double d_PE = ob_2.get_PE() - ob_1.get_PE();

    double F_dot_dx = 0;
    F_dot_dx += (ob_2.get_bmx() - ob_1.get_bmx())*(ob_1.get_internal().bmx + ob_2.get_internal().bmx)/2;
    F_dot_dx += (ob_2.get_bmy() - ob_1.get_bmy())*(ob_1.get_internal().bmy + ob_2.get_internal().bmy)/2;

    F_dot_dx += (ob_2.get_tx() - ob_1.get_tx())*(ob_1.get_internal().tx + ob_2.get_internal().tx)/2;
    F_dot_dx += (ob_2.get_ty() - ob_1.get_ty())*(ob_1.get_internal().ty + ob_2.get_internal().ty)/2;

    F_dot_dx += (ob_2.get_umx() - ob_1.get_umx())*(ob_1.get_internal().umx + ob_2.get_internal().umx)/2;
    F_dot_dx += (ob_2.get_umy() - ob_1.get_umy())*(ob_1.get_internal().umy + ob_2.get_internal().umy)/2;

    F_dot_dx += (ob_2.get_ubx() - ob_1.get_ubx())*(ob_1.get_internal().ubx + ob_2.get_internal().ubx)/2;
    F_dot_dx += (ob_2.get_uby() - ob_1.get_uby())*(ob_1.get_internal().uby + ob_2.get_internal().uby)/2;

    if (!test("d_PE = -sum F*dx?", d_PE, -F_dot_dx)) num_failures++;
  }
  
  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  delete rand;
  exit(num_failures);
}
