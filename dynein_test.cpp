#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"

#define EPSILON 1e-12
int runtime = 0;
double dt = 0.1;

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
  // Dynein in normal pentagon conformation, check if velocities agree with definitions.

  //runtime = 1*runtime;

  double bla_eq = pre_powerstroke_leftbound_internal_angles.bla;
  double la_eq  = pre_powerstroke_leftbound_internal_angles.la;
  double ta_eq  = pre_powerstroke_leftbound_internal_angles.ta;
  double ra_eq  = pre_powerstroke_leftbound_internal_angles.ra;
  
  forces no_forces =    {0,0,0,0,0,0,0,0,0,0}; // blx, bly, mlx, mly, ...
  forces right_forces = {1,0,1,0,1,0,1,0,1,0};
  forces left_forces =  {-1,0,-1,0,-1,0,-1,0,-1,0};
  forces up_forces =    {0,1,0,1,0,1,0,1,0,1};
  
  int num_failures = 0;
  {
    Dynein* dyn = new Dynein(bla_eq,                                           // starting bla
                             bla_eq + la_eq - M_PI,                            // starting mla
                             bla_eq + la_eq - M_PI + ta_eq,                    // starting mra
                             bla_eq + la_eq - M_PI + ta_eq + M_PI - ra_eq,     // starting bra
			     LEFTBOUND,                                        // starting state
			     &no_forces,                         // optional specified internal forces
			     &no_forces,                         // optional specified brownian forces
		             (equilibrium_angles*) NULL);        // optional specified equilibrium angles

    printf("Test: Dynein pentagon conformation, no internal forces, no Brownian forces.\n");
    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);
    
    num_failures += test("Is d_bla zero", dyn->get_d_bla(), 0);
    num_failures += test("Is d_mla zero", dyn->get_d_mla(), 0);
    num_failures += test("Is d_mra zero", dyn->get_d_mra(), 0);
    num_failures += test("Is d_bra zero", dyn->get_d_bra(), 0);    
    
    free(dyn);
  }

  {
    Dynein* dyn = new Dynein((90.0 / 180) * M_PI,
                             (90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
                             -(90.0 / 180) * M_PI,
			     LEFTBOUND,
			     &no_forces,
			     &no_forces,
			     NULL);
    
    printf("\nTest: Dynein vertical conformation, no internal forces, no Brownian forces.\n");
    num_failures += test("Is brx zero", dyn->get_brx(), 0);
    num_failures += test("Is mrx zero", dyn->get_mrx(), 0);
    num_failures += test("Is tx zero", dyn->get_tx(), 0);
    num_failures += test("Is mlx zero", dyn->get_mlx(), 0);
    num_failures += test("Is blx zero", dyn->get_blx(), 0);
    
    num_failures += test("Is bry fully extended", dyn->get_bry(), 2*ls + 2*lt);
    num_failures += test("Is mry fully extended", dyn->get_mry(), 1*ls + 2*lt);
    num_failures += test("Is ty fully extended", dyn->get_ty(), 1*ls + 1*lt);
    num_failures += test("Is mly fully extended", dyn->get_mly(), 1*ls);
    num_failures += test("Is bly fully extended", dyn->get_bly(), 0);
    
    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);

    num_failures += test("Is d_bla zero", dyn->get_d_bla(), 0);
    num_failures += test("Is d_mla zero", dyn->get_d_mla(), 0);
    num_failures += test("Is d_mra zero", dyn->get_d_mra(), 0);
    num_failures += test("Is d_bra zero", dyn->get_d_bra(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein((0.0 / 180) * M_PI,
                             (0.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
                             -(180.0 / 180) * M_PI,
			     LEFTBOUND,
			     &no_forces,
			     &no_forces,
			     NULL);
    
    printf("\nTest: Dynein horizontal conformation, no internal forces, no Brownian forces.\n");
    num_failures += test("Is bry zero", dyn->get_bry(), 0);
    num_failures += test("Is mry zero", dyn->get_mry(), 0);
    num_failures += test("Is ty zero", dyn->get_ty(), 0);
    num_failures += test("Is mly zero", dyn->get_mly(), 0);
    num_failures += test("Is bly zero", dyn->get_bly(), 0);

    num_failures += test("Is brx fully extended", dyn->get_brx(), 2*ls + 2*lt);
    num_failures += test("Is mrx fully extended", dyn->get_mrx(), 1*ls + 2*lt);
    num_failures += test("Is tx fully extended", dyn->get_tx(), 1*ls + 1*lt);
    num_failures += test("Is mlx fully extended", dyn->get_mlx(), 1*ls);
    num_failures += test("Is blx fully extended", dyn->get_blx(), 0);

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);

    num_failures += test("Is d_bla zero", dyn->get_d_bla(), 0);
    num_failures += test("Is d_mla zero", dyn->get_d_mla(), 0);
    num_failures += test("Is d_mra zero", dyn->get_d_mra(), 0);
    num_failures += test("Is d_bra zero", dyn->get_d_bra(), 0);
    
    free(dyn);
  }
  
  {
    Dynein* dyn = new Dynein(0.5*M_PI, 0.5*M_PI, 0.5*M_PI, 0.5*M_PI,
			     LEFTBOUND,
			     &no_forces,
			     &right_forces,
			     NULL);

    printf("\nTest: Dynein vertical, no internal forces, Brownian forces in positive x direction.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test_noteq("Is d_mlx nonzero", dyn->get_d_mlx(), 0);
    num_failures += test_noteq("Is d_tx nonzero", dyn->get_d_tx(), 0);
    num_failures += test_noteq("Is d_mrx nonzero", dyn->get_d_mrx(), 0);
    num_failures += test_noteq("Is d_brx nonzero", dyn->get_d_brx(), 0);

    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein(0.5*M_PI, 0.5*M_PI, 0.5*M_PI, 0.5*M_PI,
			     LEFTBOUND,
			     &no_forces,
			     &up_forces,
			     NULL);

    printf("\nTest: Dynein vertical, no internal forces, Brownian forces in positive y direction.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);

    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein(0, 0, 0, 0,
			     LEFTBOUND,
			     &no_forces,
			     &right_forces,
			     NULL);

    printf("\nTest: Dynein horizontal, no internal forces, Brownian forces in positive x direction.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);

    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);

    free(dyn);
  }


  {
    Dynein* dyn = new Dynein(0, 0, 0, 0,
			     LEFTBOUND,
			     &no_forces,
			     &up_forces,
			     NULL);

    printf("\nTest: Dynein horizontal, no internal forces, Brownian forces in positive y direction.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);

    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test_noteq("Is d_mly nonzero", dyn->get_d_mly(), 0);
    num_failures += test_noteq("Is d_ty nonzero", dyn->get_d_ty(), 0);
    num_failures += test_noteq("Is d_mry nonzero", dyn->get_d_mry(), 0);
    num_failures += test_noteq("Is d_bry nonzero", dyn->get_d_bry(), 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein(bla_eq,
                             bla_eq + la_eq - M_PI,
                             bla_eq + la_eq - M_PI + ta_eq,
                             bla_eq + la_eq - M_PI + ta_eq + M_PI - ra_eq,
                             LEFTBOUND,
			     NULL,
			     &no_forces,
			     NULL);

    printf("\nTest: Dynein equilibrium, checking internal forces, no Brownian forces.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_tx zero", dyn->get_d_tx(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);

    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_ty zero", dyn->get_d_ty(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);

    num_failures += test("Is f_blx zero", dyn->get_internal().blx, 0);
    num_failures += test("Is f_mlx zero", dyn->get_internal().mlx, 0);
    num_failures += test("Is f_tx zero",  dyn->get_internal().tx, 0);
    num_failures += test("Is f_mrx zero", dyn->get_internal().mrx, 0);
    num_failures += test("Is f_brx zero", dyn->get_internal().brx, 0);

    num_failures += test("Is f_bly zero", dyn->get_internal().bly, 0);
    num_failures += test("Is f_mly zero", dyn->get_internal().mly, 0);
    num_failures += test("Is f_ty zero",  dyn->get_internal().ty, 0);
    num_failures += test("Is f_mry zero", dyn->get_internal().mry, 0);
    num_failures += test("Is f_bry zero", dyn->get_internal().bry, 0);

    free(dyn);
  }

  {
    Dynein* dyn = new Dynein((108.0 / 180) * M_PI,
                             (36.0 / 180) * M_PI,
                             (144.0 / 180) * M_PI,
                             (72.0 / 180) * M_PI,
			     LEFTBOUND,
			     &no_forces,
			     &right_forces,
			     NULL);
    
    printf("\nTest: Dynein pentagon conformation, no internal forces, Brownian forces in positive x direction.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test_noteq("Is d_mlx nonzero", dyn->get_d_mlx(), 0);
    num_failures += test_noteq("Is d_mrx nonzero", dyn->get_d_mrx(), 0);
    num_failures += test_noteq("Is d_brx nonzero", dyn->get_d_brx(), 0);
    
    free(dyn);
  }
  
  {
    Dynein* dyn = new Dynein((108.0 / 180) * M_PI,
                             (36.0 / 180) * M_PI,
                             (144.0 / 180) * M_PI,
                             (72.0 / 180) * M_PI,
			     LEFTBOUND,
			     &left_forces,
			     &no_forces,
			     NULL);
    
    printf("\nTest: Dynein pentagon conformation, internal forces in negative x direction, no Brownian forces.\n");

    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test_noteq("Is d_mlx nonzero", dyn->get_d_mlx(), 0);
    num_failures += test_noteq("Is d_mrx nonzero", dyn->get_d_mrx(), 0);
    num_failures += test_noteq("Is d_brx nonzero", dyn->get_d_brx(), 0);
    
    free(dyn);
  }
  
  {
    Dynein* dyn = new Dynein((108.0 / 180) * M_PI,
                             (36.0 / 180) * M_PI,
                             (144.0 / 180) * M_PI,
                             (72.0 / 180) * M_PI,
  			     LEFTBOUND,
  			     NULL,
			     &no_forces,
			     NULL);

    printf("Test: Dynein pentagon conformation, internal forces, no Brownian forces.\n");
    num_failures += test("Is d_blx zero", dyn->get_d_blx(), 0);
    num_failures += test("Is d_bly zero", dyn->get_d_bly(), 0);
    num_failures += test("Is d_mlx zero", dyn->get_d_mlx(), 0);
    num_failures += test("Is d_mly zero", dyn->get_d_mly(), 0);
    num_failures += test("Is d_mrx zero", dyn->get_d_mrx(), 0);
    num_failures += test("Is d_mry zero", dyn->get_d_mry(), 0);
    num_failures += test("Is d_brx zero", dyn->get_d_brx(), 0);
    num_failures += test("Is d_bry zero", dyn->get_d_bry(), 0);
    
    num_failures += test("Is d_bla zero", dyn->get_d_bla(), 0);
    num_failures += test("Is d_mla zero", dyn->get_d_mla(), 0);
    num_failures += test("Is d_mra zero", dyn->get_d_mra(), 0);
    num_failures += test("Is d_bra zero", dyn->get_d_bra(), 0);    
    
    free(dyn);
  }

  
  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  exit(num_failures);
}
