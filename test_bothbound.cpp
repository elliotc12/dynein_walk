#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fenv.h>
#include <csignal>

#include "dynein_struct.h"

/*
 * Test dynein_bothbound state for errors in position, velocity and force calculations.
 */

double dt = 1e-12; // s

double kb = 1.3806e-5; // nm^2 * kg / (s^2 * K)
double T = 293.0; // K

double lt = 15.0;
double ls = 21.22;
double Lt = 15.0;  // nm, guess - not sure how DNA tail-bridge works
double Ls = 21.22; // nm, derived from PyMol dynein crystal struct 3VKH, 212.2 angstroms

// tail domain radius, not sure how to get since no info on DNA
// tail-bridge
double fake_radius_t = 1.5;  // nm
// motor domain radius, derived from PyMol, motor radius 148.6
// angstroms
double fake_radius_m = 1.48; // nm
// binding domain radius, derived from PyMol, binding radius 14.78
// angstroms
double fake_radius_b = 0.14; // nm

// water viscosity is about 0.7 mPa s = 7e-4 Pa s = 7e-4 m^2 * kg/s / m^3
// mu = 7e-4 kg/s / m
//  ... thus this is ...
// mu = 7e-4 kg/s / (1e9 nm)
//    = 7e-13 kg/s / nm
double water_viscosity_mu = 7e-13; // kg/(s*nm)

double gt = fake_radius_t*6*M_PI*water_viscosity_mu; // kg / s
double gm = fake_radius_m*6*M_PI*water_viscosity_mu; // kg / s
double gb = fake_radius_b*6*M_PI*water_viscosity_mu; // kg / s

double ct = 0.001; // force*distance = energy = nm^2 * kg / s^2
double cm = 0.5; // ???
double cb = 0.5; // ???

double ONEBOUND_UNBINDING_FORCE = 8e12;
double BOTHBOUND_UNBINDING_FORCE = 2e12;

double MICROTUBULE_REPULSION_FORCE = 30.0; // N/nm

double MICROTUBULE_BINDING_DISTANCE = 0.2; // nm

double RAND_INIT_SEED = 0;

onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.0 * M_PI,
  0.6 * M_PI
};

bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  0.5 * M_PI,
  0.5 * M_PI,
  1.0 * M_PI,
  1.5 * M_PI,
  0.5 * M_PI
};

double EPSILON = 1e-12;

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
    printf("%45s: FAIL! %g != %g.\n", msg, one, two);
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

int main(int argvc, char **argv) {
  if (FP_EXCEPTION_FATAL) {
    feenableexcept(FE_ALL_EXCEPT); // NaN generation kills program
    signal(SIGFPE, FPE_signal_handler);
  }

  double R = sqrt(2*kb*T/(gm*dt)); // Brownian force constant
  MTRand* rand = new MTRand(RAND_INIT_SEED);
  int num_failures = 0;
  
  bothbound_forces no_forces  = {0,0,0,0,0,0,0,0,0,0}; // bbx, bby, bmx, bmy, ...
  bothbound_forces out_forces = {0,0,-R,0,0,0,R,0,0,0};
  bothbound_forces x_forces   = {0,0,R,0,R,0,R,0,0,0};

  printf("****************Starting Bothbound Test****************\n\n");

  { printf("**Upwards line conformation with no forces**\n");

    bothbound_equilibrium_angles line_eq_angles = {
      M_PI/2, M_PI, 0, M_PI, M_PI/2
    };

    Dynein_bothbound dyn_bb(M_PI - 1e-15,      // nma_init
                            M_PI + 1e-15,      // fma_init
                            0,                 // nbx_init
                            0,                 // nby_init
                            1e-25,             // L
                            NULL,              // internal forces
			    &no_forces,        // brownian forces
			    &line_eq_angles,   // equilibrium angles
			    rand);             // MTRand

    printf("\tTesting coordinate system:\n");
    if (!test("nbx zero?", dyn_bb.get_nbx(), 0, 1e-6)) num_failures++;
    if (!test("nmx zero?", dyn_bb.get_nmx(), 0, 1e-6)) num_failures++;
    if (!test("tx  zero?", dyn_bb.get_tx(),  0, 1e-6)) num_failures++;
    if (!test("fmx zero?", dyn_bb.get_fmx(), 0, 1e-6)) num_failures++;
    if (!test("fbx zero?", dyn_bb.get_fbx(), 0, 1e-6)) num_failures++;

    if (!test("nby zero?", dyn_bb.get_nby(), 0, 1e-6)) num_failures++;
    if (!test("nmy = Ls?", dyn_bb.get_nmy(), Ls, 1e-6)) num_failures++;
    if (!test("ty = Ls+Lt?", dyn_bb.get_ty(), Ls+Lt, 1e-6)) num_failures++;
    if (!test("fmy = Ls?", dyn_bb.get_fmy(), Ls, 1e-6)) num_failures++;
    if (!test("fby zero?", dyn_bb.get_fby(), 0, 1e-6)) num_failures++;

    if (!test("nba - eq.nba zero?", dyn_bb.get_nba() - line_eq_angles.nba, 0, 1e-6)) num_failures++;
    if (!test("nma - eq.nma zero?", dyn_bb.get_nma() - line_eq_angles.nma, 0, 1e-6)) num_failures++;
    if (!test("ta - eq.ta zero?", dyn_bb.get_fma() + dyn_bb.get_fba() - dyn_bb.get_nma()
	      - dyn_bb.get_nba() - line_eq_angles.ta, 0, 1e-6)) num_failures++;
    if (!test("fma - eq.fma zero?", dyn_bb.get_fma() - line_eq_angles.fma, 0, 1e-6)) num_failures++;
    if (!test("fba - eq.fba zero?", dyn_bb.get_fba() - line_eq_angles.fba, 0, 1e-6)) num_failures++;

    printf("\n\tTesting force definitions:\n");
    if (!test("f.nbx zero?", dyn_bb.get_internal().nbx, 0, 1e4)) num_failures++;
    if (!test("f.nmx zero?", dyn_bb.get_internal().nmx, 0, 1e4)) num_failures++;
    if (!test("f.tx  zero?", dyn_bb.get_internal().tx , 0, 1e4)) num_failures++;
    if (!test("f.fmx zero?", dyn_bb.get_internal().fmx, 0, 1e4)) num_failures++;
    if (!test("f.fbx zero?", dyn_bb.get_internal().fbx, 0, 1e4)) num_failures++;
    if (!test("f.nby zero?", dyn_bb.get_internal().nby, 0, 1e4)) num_failures++;
    if (!test("f.nmy zero?", dyn_bb.get_internal().nmy, 0, 1e4)) num_failures++;
    if (!test("f.ty  zero?", dyn_bb.get_internal().ty , 0, 1e4)) num_failures++;
    if (!test("f.fmy zero?", dyn_bb.get_internal().fmy, 0, 1e4)) num_failures++;
    if (!test("f.fby zero?", dyn_bb.get_internal().fby, 0, 1e4)) num_failures++;
  }

   { printf("\n**House conformation with equilateral roof, outwards forces**\n");

     Dynein_bothbound dyn_bb(5*M_PI/6,     // nma_init
			     7*M_PI/6,      // fma_init
			     0,             // nbx_init
			     0,             // nby_init
			     Lt,            // L -- equilateral roof
			     &no_forces,    // internal forces
			     &out_forces,   // brownian forces
			     NULL,          // equilibrium angles
			     rand);         // MTRand                 

     printf("\tTesting coordinate definitions:\n");
     if (!test("nmx = zero?", dyn_bb.get_nmx(), 0)) num_failures++;
     if (!test("fmx = Lt?", dyn_bb.get_fmx(), Lt)) num_failures++;
     if (!test("nmy = Ls?", dyn_bb.get_nmy(), Ls)) num_failures++;
     if (!test("fmy = Ls?", dyn_bb.get_fmy(), Ls)) num_failures++;
     if (!test("tx = Lt/2?", dyn_bb.get_tx(), Lt/2)) num_failures++;
     if (!test("ty = Ls + sqrt(3)/2*Lt?", dyn_bb.get_ty(), Ls + sqrt(3)/2*Lt)) num_failures++;
     if (!test("nba = M_PI - fba?", dyn_bb.get_nba(), M_PI - dyn_bb.get_fba())) num_failures++;
     if (!test("Ln = Lf?", dyn_bb.get_ln(), dyn_bb.get_lf())) num_failures++;

     printf("\n\tTesting intermediate variable spacial derivatives:\n");
     if (!test_greater("dcosAn_dLn > 0?", dyn_bb.dcosAn_dLn, 0)) num_failures++;
     if (!test_less("dsinAn_dLn > 0?", dyn_bb.dsinAn_dLn, 0)) num_failures++;
     if (!test_less("dcosAn_dLf > 0?", dyn_bb.dcosAn_dLf, 0)) num_failures++;
     if (!test_greater("dsinAn_dLf > 0?", dyn_bb.dsinAn_dLf, 0)) num_failures++;

     if (!test_greater("dcosAns_dLn > 0?", dyn_bb.dcosAns_dLn, 0)) num_failures++;
     if (!test_less("dsinAns_dLn > 0?", dyn_bb.dsinAns_dLn, 0)) num_failures++;

     if (!test("dcosAn_dLn = -dcosAf_dLf?", dyn_bb.dcosAn_dLn, -dyn_bb.dcosAf_dLf)) num_failures++;
     if (!test("dcosAn_dLf = -dcosAf_dLn?", dyn_bb.dcosAn_dLf, -dyn_bb.dcosAf_dLn)) num_failures++;
     if (!test("dcosAns_dLn = dcosAfs_dLf?", dyn_bb.dcosAns_dLn, dyn_bb.dcosAfs_dLf)) num_failures++;
     if (!test("dcosAns_dLf = dcosAfs_dLn?", dyn_bb.dcosAns_dLf, dyn_bb.dcosAfs_dLn)) num_failures++;

     if (!test("dsinAn_dLn = dsinAf_dLf?", dyn_bb.dsinAn_dLn, dyn_bb.dsinAf_dLf)) num_failures++;
     if (!test("dsinAn_dLf = dsinAf_dLn?", dyn_bb.dsinAn_dLf, dyn_bb.dsinAf_dLn)) num_failures++;
     if (!test("dsinAns_dLn = -dsinAfs_dLf?", dyn_bb.dsinAns_dLn, -dyn_bb.dsinAfs_dLf)) num_failures++;
     if (!test("dsinAns_dLf = -dsinAfs_dLn?", dyn_bb.dsinAns_dLf, -dyn_bb.dsinAfs_dLn)) num_failures++;

     printf("\n\tTesting cartesian spacial derivatives:\n");
     if (!test("dXnm_dLn = -dXfm_dLf?", dyn_bb.dXnm_dLn, -dyn_bb.dXfm_dLf)) num_failures++;
     if (!test("dXnm_dLf = -dXfm_dLn?", dyn_bb.dXnm_dLf, -dyn_bb.dXfm_dLn)) num_failures++;
     if (!test("dYnm_dLn = dYfm_dLf?", dyn_bb.dYnm_dLn, dyn_bb.dYfm_dLf)) num_failures++;
     if (!test("dYnm_dLf = dYfm_dLn?", dyn_bb.dYnm_dLf, dyn_bb.dYfm_dLn)) num_failures++;

     if (!test("dXt_dLn = -dXt_dLf?", dyn_bb.dXt_dLn, -dyn_bb.dXt_dLf)) num_failures++;
     if (!test("dYt_dLn = dYt_dLf?", dyn_bb.dYt_dLn, dyn_bb.dYt_dLf)) num_failures++;

     printf("\n\tTesting time derivatives:\n");
     if (!test("d_Ln = d_Lf?", dyn_bb.get_d_Ln(), dyn_bb.get_d_Lf())) num_failures++;
     if (!test("d_nma = -d_fma?", dyn_bb.get_d_nma(), -dyn_bb.get_d_fma())) num_failures++;
     if (!test_less("d_nmx < 0?", dyn_bb.get_d_nmx(), 0)) num_failures++;
     if (!test("d_tx = zero?", dyn_bb.get_d_tx(), 0, 1e-5)) num_failures++;
     if (!test_less("d_ty < 0?", dyn_bb.get_d_ty(), 0)) num_failures++;
     if (!test_greater("d_fmx > 0?", dyn_bb.get_d_fmx(), 0)) num_failures++;
     if (!test("d_nmx = -d_fmx?", dyn_bb.get_d_nmx(), -dyn_bb.get_d_fmx())) num_failures++;
     if (!test("d_nmy = d_fmy?", dyn_bb.get_d_nmy(), dyn_bb.get_d_fmy())) num_failures++;
  }

  { printf("\n**'Almost' upwards line conformation with +x forces**\n");

    Dynein_bothbound dyn_bb(M_PI - 1e-7,      // nma_init
                            M_PI + 1e-7,      // fma_init
                            0,                // nbx_init
                            0,                // nby_init
                            1e-1,             // L
                            &no_forces,       // internal forces
			    &x_forces,        // brownian forces
			    NULL,             // equilibrium angles
			    rand);            // MTRand            

    printf("\tTesting time derivatives:\n");
    if (!test_noteq("d_nmx_dt nonzero?", dyn_bb.get_d_nmx(), 0)) num_failures++;
    if (!test_noteq("d_tx_dt  nonzero?", dyn_bb.get_d_tx(), 0)) num_failures++;
    if (!test_noteq("d_fmx_dt nonzero?", dyn_bb.get_d_fmx(), 0)) num_failures++;
  }

  { printf("\n**Two table-ish conformations with near/far domains flipped**\n");

    double t_nma = acos(Lt/(2*Ls));
    
    bothbound_equilibrium_angles left_table_eq_angles = {
      M_PI - 2*t_nma + M_PI/3,
      t_nma,
      2*t_nma + M_PI/3,
      2*M_PI - t_nma,
      2*t_nma - M_PI/3
    };

    bothbound_equilibrium_angles right_table_eq_angles = {
      2*t_nma - M_PI/3,
      2*M_PI - t_nma,
      -(2*t_nma + M_PI/3),
      t_nma,
      M_PI - 2*t_nma + M_PI/3
    };

    Dynein_bothbound left_dyn_bb(t_nma,                  // nma_init	       	
				 2*M_PI - t_nma,         // fma_init	       	
				 0,                      // nbx_init	       	
				 0,                      // nby_init	       	
				 Ls,                     // L		       	
				 NULL,                   // internal forces    
				 &no_forces,             // brownian forces    
				 &left_table_eq_angles,  // equilibrium angles 
				 rand);                  // MTRand             

    Dynein_bothbound right_dyn_bb(2*M_PI - t_nma,         // nma_init	       	
				  t_nma,                  // fma_init	       	
				  Ls,                     // nbx_init	       	
				  0,                      // nby_init	       	
				  -Ls,                    // L		       	
				  NULL,                   // internal forces    
				  &no_forces,             // brownian forces    
				  &right_table_eq_angles, // equilibrium angles 
				  rand);                  // MTRand             

    printf("\tTesting coordinate definitions:\n");
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

    if (!test("left by coords equal?",
	      left_dyn_bb.get_nby(), right_dyn_bb.get_fby())) num_failures++;
    if (!test("left my coords equal?",
	      left_dyn_bb.get_nmy(), right_dyn_bb.get_fmy())) num_failures++;
    if (!test("ty coords equal?",
	      left_dyn_bb.get_ty(), right_dyn_bb.get_ty())) num_failures++;
    if (!test("right my coords equal?",
	      left_dyn_bb.get_fmy(), right_dyn_bb.get_nmy())) num_failures++;
    if (!test("right by coords equal?",
	      left_dyn_bb.get_fby(), right_dyn_bb.get_nby())) num_failures++;
    
    if (!test("Left Ln/Lf equal?",
	      left_dyn_bb.get_ln(), right_dyn_bb.get_lf())) num_failures++;
    if (!test("Right Ln/Lf equal?",
	      left_dyn_bb.get_lf(), right_dyn_bb.get_ln())) num_failures++;

    if (!test("Left binding angle equal?",
	      left_dyn_bb.get_nba(), right_dyn_bb.get_fba())) num_failures++;
    if (!test("Right binding angle equal?",
	      left_dyn_bb.get_fba(), right_dyn_bb.get_nba())) num_failures++;

    printf("\n\tTesting internal forces:\n");
    if (!test("left f.nmx = right f.fmx?", left_dyn_bb.get_internal().nmx,
	      right_dyn_bb.get_internal().fmx)) num_failures++;
    if (!test("left f.tx = right f.tx?", left_dyn_bb.get_internal().tx,
	      right_dyn_bb.get_internal().tx)) num_failures++;
    if (!test("left f.fmx = right f.nmx?", left_dyn_bb.get_internal().fmx,
	      right_dyn_bb.get_internal().nmx)) num_failures++;
    if (!test("left f.nmy = right f.fmy?", left_dyn_bb.get_internal().nmy,
	      right_dyn_bb.get_internal().fmy)) num_failures++;
    if (!test("left f.ty = right f.ty?", left_dyn_bb.get_internal().ty,
	      right_dyn_bb.get_internal().ty)) num_failures++;
    if (!test("left f.fmy = right f.nmy?", left_dyn_bb.get_internal().fmy,
	      right_dyn_bb.get_internal().nmy)) num_failures++;
    
    printf("\n\tTesting time derivatives:\n");
    if (!test("left angle velocities exactly equal?",
	      left_dyn_bb.get_d_nma(), right_dyn_bb.get_d_fma(), 0)) num_failures++;
    if (!test("right angle velocities exactly equal?",
	      left_dyn_bb.get_d_fma(), right_dyn_bb.get_d_nma(), 0)) num_failures++;

    if (!test("left mx velocities exactly equal?",
	      left_dyn_bb.get_d_nmx(), right_dyn_bb.get_d_fmx(), 0)) num_failures++;
    if (!test("tx velocities exactly equal?",
	      left_dyn_bb.get_d_tx(), right_dyn_bb.get_d_tx(), 0)) num_failures++;
    if (!test("right mx velocities exactly equal?",
	      left_dyn_bb.get_d_fmx(), right_dyn_bb.get_d_nmx(), 0)) num_failures++;

    if (!test("left my velocities exactly equal?",
	      left_dyn_bb.get_d_nmy(), right_dyn_bb.get_d_fmy(), 0)) num_failures++;
    if (!test("ty velocities exactly equal?",
	      left_dyn_bb.get_d_ty(), right_dyn_bb.get_d_ty(), 0)) num_failures++;
    if (!test("right my velocities exactly equal?",
	      left_dyn_bb.get_d_fmy(), right_dyn_bb.get_d_nmy(), 0)) num_failures++;
  }

  {

    double poke = 0.01;
    Dynein_bothbound bb_1(5*M_PI/6,      // nma_init
			  7*M_PI/6,      // fma_init
			  0,             // nbx_init
			  0,             // nby_init
			  Lt,            // L -- equilateral roof
			  NULL,          // internal forces
			  &no_forces,    // brownian forces
			  NULL,          // equilibrium angles
			  rand);         // MTRand                 

    Dynein_bothbound bb_2(5*M_PI/6 + poke, // nma_init
			  7*M_PI/6 + poke, // fma_init
			  0,               // nbx_init
			  0,               // nby_init
			  Lt,              // L -- equilateral roof
			  NULL,            // internal forces
			  &no_forces,      // brownian forces
			  NULL,            // equilibrium angles
			  rand);           // MTRand                 

    double d_PE = bb_2.get_PE() - bb_1.get_PE();

    double d_nb = sqrt((bb_2.get_nbx() - bb_1.get_nbx())*(bb_2.get_nbx() - bb_1.get_nbx())
		       + (bb_2.get_nby() - bb_1.get_nby())*(bb_2.get_nby() - bb_1.get_nby()));
    double d_nm = sqrt((bb_2.get_nmx() - bb_1.get_nmx())*(bb_2.get_nmx() - bb_1.get_nmx())
		       + (bb_2.get_nmy() - bb_1.get_nmy())*(bb_2.get_nmy() - bb_1.get_nmy()));
    double d_t = sqrt((bb_2.get_tx() - bb_1.get_tx())*(bb_2.get_tx() - bb_1.get_tx())
		       + (bb_2.get_ty() - bb_1.get_ty())*(bb_2.get_ty() - bb_1.get_ty()));
    double d_fm = sqrt((bb_2.get_fmx() - bb_1.get_fmx())*(bb_2.get_fmx() - bb_1.get_fmx())
		       + (bb_2.get_fmy() - bb_1.get_fmy())*(bb_2.get_fmy() - bb_1.get_fmy()));
    double d_fb = sqrt((bb_2.get_fbx() - bb_1.get_fbx())*(bb_2.get_fbx() - bb_1.get_fbx())
		       + (bb_2.get_fby() - bb_1.get_fby())*(bb_2.get_fby() - bb_1.get_fby()));

    double f_nb = sqrt(bb_1.get_internal().nbx*bb_1.get_internal().nbx
		       + bb_1.get_internal().nby*bb_1.get_internal().nby);
    double f_nm = sqrt(bb_1.get_internal().nmx*bb_1.get_internal().nmx
		       + bb_1.get_internal().nmy*bb_1.get_internal().nmy);
    double f_t = sqrt(bb_1.get_internal().tx*bb_1.get_internal().tx
		       + bb_1.get_internal().ty*bb_1.get_internal().ty);
    double f_fm = sqrt(bb_1.get_internal().fmx*bb_1.get_internal().fmx
		       + bb_1.get_internal().fmy*bb_1.get_internal().fmy);
    double f_fb = sqrt(bb_1.get_internal().fbx*bb_1.get_internal().fbx
		       + bb_1.get_internal().fby*bb_1.get_internal().fby);

    if (!test("d_PE = sum F*dx?",
	      d_PE, d_nb*f_nb + d_nm*f_nm + d_t*f_t + d_fm*f_fm + d_fb*f_fb)) num_failures++;
  }

  if (num_failures == 0) {
    printf("All %d tests pass!\n\n", num_tests);
  } else {
    printf("%d/%d tests fail!\n\n", num_failures, num_tests);
  }
  exit(num_failures);
}
