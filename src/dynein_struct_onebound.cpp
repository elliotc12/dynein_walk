#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"

// static bool OB_PHYSICAL = true;

/* ********************* ONEBOUND DYNEIN FUNCTIONS ************************** */

Dynein_onebound::Dynein_onebound(double bba_init, double bma_init,
				 double uma_init, double uba_init,
				 double bbx_init, double bby_init,
				 State s, onebound_forces *internal_test,
				 onebound_forces *brownian_test,
				 onebound_equilibrium_angles* eq_angles,
				 Rand* mtrand) {
  bbx = bbx_init;
  bby = bby_init;

  bba = bba_init;
  bma = bma_init;
  uma = uma_init;
  uba = uba_init;

  state = s;
  internal_testcase = internal_test;
  brownian_testcase = brownian_test;

  if (eq_angles) {
    eq = *eq_angles; // use test angles
  } else {
    eq = onebound_post_powerstroke_internal_angles; // use experimental angles
  }

  rand = mtrand;

  update_velocities();
}

Dynein_onebound::Dynein_onebound(Dynein_bothbound* old_dynein, Rand* mtrand, State s) {
  if (s == NEARBOUND) {
    bbx = old_dynein->get_nbx();
    bby = 0;

    state = NEARBOUND;

    bba = old_dynein->get_nba();
    bma = old_dynein->get_nma() + old_dynein->get_nba() - M_PI;
    uma = old_dynein->get_fma() + old_dynein->get_fba() - M_PI;
    uba = old_dynein->get_fba();

  } else {
    bbx = old_dynein->get_fbx();
    bby = 0;

    state = FARBOUND;

    bba = old_dynein->get_fba();
    bma = old_dynein->get_fma() + old_dynein->get_fba() - M_PI;
    uma = old_dynein->get_nma() + old_dynein->get_nba() - M_PI;
    uba = old_dynein->get_nba();
  }

  while (bba < 0)      bba += 2*M_PI;
  while (bba > 2*M_PI) bba -= 2*M_PI;

  while (bma < 0)      bma += 2*M_PI;
  while (bma > 2*M_PI) bma -= 2*M_PI;

  while (uma < 0)      uma += 2*M_PI;
  while (uma > 2*M_PI) uma -= 2*M_PI;

  while (uba < 0)      uba += 2*M_PI;
  while (uba > 2*M_PI) uba -= 2*M_PI;

  internal_testcase = NULL;
  brownian_testcase = NULL;

  eq = onebound_post_powerstroke_internal_angles; // use experimental angles

  rand = mtrand;

  update_velocities();

  if (am_debugging_conversions) {
    printf("DEBUG:\nDEBUG: creating onebound from bothbound!\n");
    if (get_state() == NEARBOUND) {
      printf("DEBUG: nbx/bbx = %8g vs %8g  nby/bby = %8g vs %8g\n",
	     old_dynein->get_nbx(), get_bbx(), old_dynein->get_nby(), get_bby());
      printf("DEBUG: nmx/bmx = %8g vs %8g  nmy/bmy = %8g vs %8g\n",
	     old_dynein->get_nmx(), get_bmx(), old_dynein->get_nmy(), get_bmy());
      printf("DEBUG: fbx/ubx = %8g vs %8g  fby/uby = %8g vs %8g\n",
	     old_dynein->get_fbx(), get_ubx(), old_dynein->get_fby(), get_uby());
    } else {
      printf("DEBUG: nbx/ubx = %8g vs %8g  nby/uby = %8g vs %8g\n",
	     old_dynein->get_nbx(), get_ubx(), old_dynein->get_nby(), get_uby());
      printf("DEBUG: nmx/umx = %8g vs %8g  nmy/umy = %8g vs %8g\n",
	     old_dynein->get_nmx(), get_umx(), old_dynein->get_nmy(), get_umy());
      printf("DEBUG: fbx/bbx = %8g vs %8g  fby/bby = %8g vs %8g\n",
	     old_dynein->get_fbx(), get_bbx(), old_dynein->get_fby(), get_bby());
    }
    printf("DEBUG:      tx = %8g vs %8g       ty = %8g vs %8g\n",
	   old_dynein->get_tx(), get_tx(), old_dynein->get_ty(), get_ty());
  }
}

void Dynein_onebound::update_brownian_forces() {
  if (brownian_testcase) {
    r = *brownian_testcase; // just copy over forces!
  } else {

    rand->gauss2(sqrt(2*kb*T*gb/dt), &r.bbx, &r.bby);
    rand->gauss2(sqrt(2*kb*T*gm/dt), &r.bmx, &r.bmy);
    rand->gauss2(sqrt(2*kb*T*gt/dt), &r.tx, &r.ty);
    rand->gauss2(sqrt(2*kb*T*gm/dt), &r.umx, &r.umy);
    rand->gauss2(sqrt(2*kb*T*gb/dt), &r.ubx, &r.uby);
  }
}

void Dynein_onebound::update_internal_forces() {
  if (internal_testcase) {
    f = *internal_testcase;
  } else {
    f.bbx = 0;     f.bby = 0;     // Initialize forces to zero
    f.bmx = 0;     f.bmy = 0;
    f.tx  = 0;     f.ty  = 0;
    f.umx = 0;     f.umy = 0;
    f.ubx = 0;     f.uby = 0;

    double T, f1, f2, f1x, f1y, f2x, f2y;

    T = cb*(bba - eq.bba);
    PE_bba = 0.5*cb*(bba - eq.bba)*(bba - eq.bba);
    f2 = T/Ls;
    f2x = f2 * sin(bba);
    f2y = f2 * -cos(bba);
    f.bmx += f2x;
    f.bmy += f2y;
    f.bbx += -f2x; // Equal and opposite forces!  :)
    f.bby += -f2y; // Equal and opposite forces!  :)

    T = cm*((bma + M_PI - bba) - eq.bma);
    PE_bma = 0.5*cm*((bma + M_PI - bba) - eq.bma)*((bma + M_PI - bba) - eq.bma);
    f1 = T/Ls;
    f2 = T/Lt;
    f1x = f1 * sin(bba);
    f1y = f1 * -cos(bba);
    f2x = f2 * sin(bma);
    f2y = f2 * -cos(bma);
    f.bbx += f1x;
    f.bby += f1y;
    f.tx  += f2x;
    f.ty  += f2y;
    f.bmx += -(f1x + f2x);
    f.bmy += -(f1y + f2y);

    T = ct*(uma - bma - eq.ta);
    PE_ta = 0.5*ct*(uma - bma - eq.ta)*(uma - bma - eq.ta);
    f1 = T / Lt;
    f2 = T / Lt;
    f1x = f1 * sin(bma);
    f1y = f1 * -cos(bma);
    f2x = f2 * -sin(uma);
    f2y = f2 * cos(uma);
    f.bmx += f1x;
    f.bmy += f1y;
    f.umx += f2x;
    f.umy += f2y;
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = cm*((uma + M_PI - uba) - eq.uma);
    PE_uma = 0.5*cm*((uma + M_PI - uba) - eq.uma)*((uma + M_PI - uba) - eq.uma);
    f1 = T / Lt;
    f2 = T / Ls;
    f1x = f1 * sin(uma);
    f1y = f1 * -cos(uma);
    f2x = f2 * sin(uba);
    f2y = f2 * -cos(uba);
    f.tx  += f1x;
    f.ty  += f1y;
    f.ubx += f2x;
    f.uby += f2y;
    f.umx += -(f1x + f2x);
    f.umy += -(f1y + f2y);

    f.tx += tail_force;

//    if (get_bmy() < 0) f.bmy += MICROTUBULE_REPULSION_FORCE * fabs(get_bmy());
//    if (get_ty()  < 0) f.ty  += MICROTUBULE_REPULSION_FORCE * fabs(get_ty());
//    if (get_umy() < 0) f.umy += MICROTUBULE_REPULSION_FORCE * fabs(get_umy());
//    if (get_uby() < 0) f.uby += MICROTUBULE_REPULSION_FORCE * fabs(get_uby());
  }
}


// void Dynein_onebound::update_velocities() {
//   update_internal_forces();
//   update_brownian_forces();

//   double A1, A2, A3, A4;  // Start solving for velocities with matrix solution in derivation.pdf
//   double B1, B2, B3, B4;
//   double C1, C2, C3, C4;
//   double D1, D2, D3, D4;
//   double X1, X2, X3, X4;
//   double Nbb, Nml, Nmr, Nbr;
//   double D;

//   A1 = -4*Ls;
//   A2 = -3*Lt*(sin(bma)*sin(bba) + cos(bma)*cos(bba));
//   A3 = 2*Lt*(sin(uma)*sin(bba) + cos(uma)*cos(bba));
//   A4 = Ls*(sin(uba)*sin(bba) + cos(uba)*cos(bba));
//   B1 = -3*Ls*(sin(bba)*sin(bma) + cos(bba)*cos(bma));
//   B2 = -3*Lt;
//   B3 = +2*Lt*(sin(uma)*sin(bma) + cos(uma)*cos(bma));
//   B4 = +Ls*(sin(uba)*sin(bma) + cos(uba)*cos(bma));
//   C1 = -2*Ls*(sin(bba)*sin(uma) + cos(bba)*cos(uma));
//   C2 = -2*Lt*(sin(bma)*sin(uma) + cos(bma)*cos(uma));
//   C3 = 2*Lt;
//   C4 = Ls*(sin(uba)*sin(uma) + cos(uba)*cos(uma));
//   D1 = Ls*(cos(uba)*cos(bba) + sin(uba)*sin(bba));
//   D2 = Lt*(cos(uba)*cos(bma) + sin(uba)*sin(bma));
//   D3 = -Lt*(cos(uba)*cos(uma) + sin(uba)*sin(uma));
//   D4 = -Ls;

//   X1 = (- 1/gm*f.bmy - 1/gt*f.ty - 1/gm*f.umy - 1/gb*f.uby - r.bmy/gm - r.ty/gt - r.umy/gm -
//	r.uby/gb)*cos(bba) +( 1/gm*f.bmx + 1/gt*f.tx + 1/gm*f.umx + 1/gb*f.ubx + r.bmx/gm +
//	r.tx/gt + r.umx/gm + r.ubx/gb )*sin(bba);

//   X2 = (-1/gt*f.ty - 1/gm*f.umy - 1/gb*f.uby - r.ty/gt - r.umy/gm - r.uby/gb)*cos(bma) +
//     (1/gt*f.tx + 1/gm*f.umx + 1/gb*f.ubx + r.tx/gt + r.umx/gm + r.ubx/gb)*sin(bma);

//   X3 = (-r.umy/gm -r.uby/gb - 1/gm*f.umy - 1/gb*f.uby)*cos(uma) + (r.umx/gm + r.ubx/gb +
//	 1/gm*f.umx + 1/gb*f.ubx)*sin(uma);

//   X4 = (r.uby/gb + 1/gb*f.uby)*cos(uba) - (r.ubx/gb + 1/gb*f.ubx)*sin(uba);

//   Nbb = (-B2*C4*D3*X1 + B2*C3*D4*X1 + A4*C3*D2*X2 - A3*C4*D2*X2 - A4*C2*D3*X2 + A2*C4*D3*X2
//	 +A3*C2*D4*X2 - A2*C3*D4*X2 + A4*B2*D3*X3 - A3*B2*D4*X3 - A4*B2*C3*X4 + A3*B2*C4*X4
//	 +B4*(-C3*D2*X1 + C2*D3*X1 + A3*D2*X3 - A2*D3*X3 - A3*C2*X4 + A2*C3*X4)
//	 +B3*(C4*D2*X1 - C2*D4*X1 - A4*D2*X3 + A2*D4*X3 + A4*C2*X4 - A2*C4*X4));

//   Nml = (B1*C4*D3*X1 - B1*C3*D4*X1 - A4*C3*D1*X2 + A3*C4*D1*X2 + A4*C1*D3*X2 - A1*C4*D3*X2
//	 -A3*C1*D4*X2 + A1*C3*D4*X2 - A4*B1*D3*X3 + A3*B1*D4*X3 + A4*B1*C3*X4 - A3*B1*C4*X4
//	 +B4*(C3*D1*X1 - C1*D3*X1 - A3*D1*X3 + A1*D3*X3 + A3*C1*X4 - A1*C3*X4)
//	 +B3*(-C4*D1*X1 + C1*D4*X1 + A4*D1*X3 - A1*D4*X3 - A4*C1*X4 + A1*C4*X4));

//   Nmr = (-B1*C4*D2*X1 + B1*C2*D4*X1 + A4*C2*D1*X2 - A2*C4*D1*X2 - A4*C1*D2*X2 + A1*C4*D2*X2
//	 +A2*C1*D4*X2 - A1*C2*D4*X2 + A4*B1*D2*X3 - A2*B1*D4*X3 - A4*B1*C2*X4 + A2*B1*C4*X4
//	 +B4*(-C2*D1*X1 + C1*D2*X1 + A2*D1*X3 - A1*D2*X3 - A2*C1*X4 + A1*C2*X4)
//	 +B2*(C4*D1*X1 - C1*D4*X1 - A4*D1*X3 + A1*D4*X3 + A4*C1*X4 - A1*C4*X4));

//   Nbr = (B1*C3*D2*X1 - B1*C2*D3*X1 - A3*C2*D1*X2 + A2*C3*D1*X2 + A3*C1*D2*X2 - A1*C3*D2*X2
//	 -A2*C1*D3*X2 + A1*C2*D3*X2 - A3*B1*D2*X3 + A2*B1*D3*X3 + A3*B1*C2*X4 - A2*B1*C3*X4
//	 +B3*(C2*D1*X1 - C1*D2*X1 - A2*D1*X3 + A1*D2*X3 + A2*C1*X4 - A1*C2*X4)
//	 +B2*(-C3*D1*X1 + C1*D3*X1 + A3*D1*X3 - A1*D3*X3 - A3*C1*X4 + A1*C3*X4));

//   D = A2*B4*C3*D1 - A2*B3*C4*D1 - A1*B4*C3*D2 + A1*B3*C4*D2 - A2*B4*C1*D3 + A1*B4*C2*D3
//          +A2*B1*C4*D3 - A1*B2*C4*D3 + A4*(B3*C2*D1 - B2*C3*D1 - B3*C1*D2 + B1*C3*D2 + B2*C1*D3
//          -B1*C2*D3)+ A2*B3*C1*D4 - A1*B3*C2*D4 - A2*B1*C3*D4 + A1*B2*C3*D4
//          +A3*(-B4*C2*D1 + B2*C4*D1 + B4*C1*D2 - B1*C4*D2 - B2*C1*D4 + B1*C2*D4);

//   //assert(D != 0);
//   if (D == 0 and OB_PHYSICAL) {
//     printf("Simulation unphysical, D == 0\n");
//     OB_PHYSICAL = false;
//   }

//   d_bba = Nbb/D;
//   d_bma = Nml/D;
//   d_uma = Nmr/D;
//   d_uba = Nbr/D;
// }

double Power(double num, int pow) {
  if (pow == 1) return num;
  else if (pow == 2) return num*num;
  else if (pow == 3) return num*num*num;
  else if (pow == 4) return (num*num)*(num*num);
  else if (pow == 5) return (num*num)*(num*num)*num;
  else if (pow == 6) return (num*num)*(num*num)*(num*num);
  else if (pow == 7) return (num*num*num)*(num*num*num)*num;
  else {
    printf("Need more power!\n");
    exit(pow);
  }
}

bool Dynein_onebound::update_velocities() {
  if (bma < -2*M_PI or bma > 2*M_PI) { // check motor angles for crazy states
    if (am_naively_correcting_nan_errors) {
      if (bma < -2*M_PI) {
        if (am_debugging_naive_corrections) printf("Naive correction: bma < -2*M_PI\n");
        bma = -2*M_PI + 1e-6;
      }
      if (bma > 2*M_PI) {
        if (am_debugging_naive_corrections) printf("Naive correction: bma > 2*M_PI\n");
        bma = 2*M_PI - 1e-6;
      }
    }
    if (am_debugging_angles) printf("bma angle is crazy man! %g\n", bma);
  }
  else if (am_debugging_angles) printf("bma angle is cool:      %g\n", bma);


  if (uma < -2*M_PI or uma > 2*M_PI) {
    if (am_naively_correcting_nan_errors) {
      if (uma < -2*M_PI) {
        if (am_debugging_naive_corrections) printf("Naive correction: uma < -2*M_PI\n");
        uma = -2*M_PI + 1e-6;
      }
      if (uma > 2*M_PI) {
        if (am_debugging_naive_corrections) printf("Naive correction: uma > 2*M_PI\n");
        uma = 2*M_PI - 1e-6;
      }
    }
    if (am_debugging_angles) printf("uma angle is crazy man! %g\n", uma);
  }
  else if (am_debugging_angles) printf("uma angle is cool:      %g\n", uma);

  if (bba > M_PI or bba < 0) { // check binding angles for crazy states
    if (am_naively_correcting_nan_errors) {
      if (bba < 0) {
        if (am_debugging_naive_corrections) printf("Naive correction: bba < 0\n");
        bba = 1e-6;
      }
      if (bba > M_PI) {
        if (am_debugging_naive_corrections) printf("Naive correction: bba > M_PI\n");
        bba = M_PI - 1e-6;
      }
    }
    if (am_debugging_angles) printf("bba angle is crazy man! %g\n", bba);
  }

  if (uba > 2*M_PI and am_naively_correcting_nan_errors) { // cyclic unbinding angle
    if (am_debugging_naive_corrections) printf("Naive correction: uba > 2*M_PI\n");
    uba = fmod(uba, 2*M_PI);
  }

  if (uba < 0 and am_naively_correcting_nan_errors) {
    if (am_debugging_naive_corrections) printf("Naive correction: uba < 0\n");
    uba = 2*M_PI - fmod(fabs(uba), 2*M_PI);
  }

  //******* Checking for sub-MT dynein ********
  if (am_crashing_on_unphysical_behavior) {
    if (get_bmy() < 0.0 or get_ty() < 0.0 or get_umy() < 0.0) {
      // printf("A domain is under the MT! bmy, ty, umy, uby: : %g, %g, %g, %g\n", get_bmy(), get_ty(), get_umy(), get_uby());
      // fprintf(stderr, "A domain is under the MT! bmy, ty, umy, uby: : %g, %g, %g, %g\n", get_bmy(), get_ty(), get_umy(), get_uby());
      // fprintf(stderr, "These are bad parameters; retrying.\n");
      if (am_only_writing_on_crash) on_crash_write_movie_buffer();
      return false;
	// exit(1);
    }
  }

  update_internal_forces();
  update_brownian_forces();

  double X_bs = get_bmx() - get_bbx();
  double X_bt = get_tx() - get_bmx();
  double X_ut = get_umx() - get_tx();
  double X_us = get_ubx() - get_umx();
  double Y_bs = get_bmy() - get_bby();
  double Y_bt = get_ty() - get_bmy();
  double Y_ut = get_umy() - get_ty();
  double Y_us = get_uby() - get_umy();

  double AA = Ls * sin(bba);
  double BB = Lt * sin(bma);
  double CC = -Lt * sin(uma);
  double DD = -Ls * sin(uba);
  double EE = -X_bs / gm;
  double FF = X_bt / gm;
  double GG = -X_bt / gt;
  double HH = X_ut / gt;
  double II = -X_ut / gm;
  double JJ = X_us / gm;
  double KK = -X_us / gb;
  double LL = -Ls * cos(bba);
  double MM = -Lt * cos(bma);
  double NN = Lt * cos(uma);
  double OO = Ls * cos(uba);
  double PP = -Y_bs / gm;
  double QQ = Y_bt / gm;
  double RR = -Y_bt / gt;
  double SS = Y_ut / gt;
  double TT = -Y_ut / gm;
  double UU = Y_us / gm;
  double VV = -Y_us / gb;

  double X1 = -(f.bmx + r.bmx) / gm;
  double X2 = -(f.tx + r.tx) / gt;
  double X3 = -(f.umx + r.umx) / gm;
  double X4 = -(f.ubx + r.ubx) / gb;
  double X5 = -(f.bmy + r.bmy) / gm;
  double X6 = -(f.ty + r.ty) / gt;
  double X7 = -(f.umy + r.umy) / gm;
  double X8 = -(f.uby + r.uby) / gb;

  d_bba = (-(GG*II*MM*NN*PP*X1) + CC*HH*MM*PP*RR*X1 -
      BB*HH*NN*PP*RR*X1 + BB*II*NN*PP*RR*X1 - CC*GG*MM*PP*SS*X1 +
      BB*GG*NN*PP*SS*X1 + CC*GG*MM*PP*TT*X1 - BB*CC*PP*RR*TT*X1 +
      FF*II*MM*NN*PP*X2 - EE*II*MM*NN*QQ*X2 + CC*FF*MM*PP*SS*X2 -
      BB*FF*NN*PP*SS*X2 - CC*EE*MM*QQ*SS*X2 + BB*EE*NN*QQ*SS*X2 -
      CC*FF*MM*PP*TT*X2 + CC*EE*MM*QQ*TT*X2 - FF*HH*MM*NN*PP*X3 +
      EE*HH*MM*NN*QQ*X3 + BB*FF*NN*PP*SS*X3 - BB*EE*NN*QQ*SS*X3 +
      EE*GG*II*MM*NN*X5 - CC*EE*HH*MM*RR*X5 + BB*EE*HH*NN*RR*X5 -
      BB*EE*II*NN*RR*X5 + CC*EE*GG*MM*SS*X5 - BB*EE*GG*NN*SS*X5 -
      CC*EE*GG*MM*TT*X5 + BB*CC*EE*RR*TT*X5 - CC*FF*HH*MM*PP*X6 +
      BB*FF*HH*NN*PP*X6 - BB*FF*II*NN*PP*X6 + CC*EE*HH*MM*QQ*X6 -
      BB*EE*HH*NN*QQ*X6 + BB*EE*II*NN*QQ*X6 + BB*CC*FF*PP*TT*X6 -
      BB*CC*EE*QQ*TT*X6 + CC*FF*HH*MM*PP*X7 - CC*EE*HH*MM*QQ*X7 -
      BB*CC*FF*PP*SS*X7 + BB*CC*EE*QQ*SS*X7)/ (EE*GG*II*LL*MM*NN +
      BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP - AA*FF*HH*MM*NN*PP +
      AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP - BB*EE*HH*LL*NN*QQ +
      BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ - AA*EE*II*MM*NN*QQ -
      CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR - BB*EE*II*LL*NN*RR +
      AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR + AA*BB*II*NN*PP*RR +
      CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS - BB*CC*FF*LL*PP*SS +
      AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS + AA*BB*GG*NN*PP*SS +
      BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS - CC*EE*GG*LL*MM*TT +
      BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT + AA*CC*GG*MM*PP*TT -
      BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT + BB*CC*EE*LL*RR*TT -
      AA*BB*CC*PP*RR*TT) - ((FF*PP - EE*QQ)*(HH*MM - BB*SS)*(JJ*NN -
      CC*UU)*(-(((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
      Power(AA,6)*Power(BB,3)*DD*HH*NN +
      Power(AA,6)*Power(BB,3)*DD*II*NN -
      Power(AA,6)*Power(BB,3)*CC*II*OO)* (-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR))) - (-(EE*LL) +
      AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
      Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
      Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
      (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS))*
      (-(((Power(AA,3)*Power(BB,2)*CC*FF*LL - Power(AA,4)*BB*CC*FF*MM
      + Power(AA,4)*BB*CC*GG*MM -
      Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
      (Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ))*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
      AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
      Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
      Power(AA,4)*Power(BB,2)*CC*X6))) +
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
      AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
      Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
      Power(AA,4)*Power(BB,2)*NN*X2 + Power(AA,4)*Power(BB,2)*NN*X3 -
      Power(AA,4)*Power(BB,2)*CC*X7)))) + (-((-(EE*LL) +
      AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
      Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
      Power(AA,4)*Power(BB,2)*GG*NN)* (-(EE*LL) + AA*PP) -
      (Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS)) + (-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN -
      Power(AA,4)*Power(BB,2)*CC*TT))*
      (-(((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
      Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
      Power(AA,6)*Power(BB,3)*DD*GG*NN)*(-(EE*LL) + AA*PP) -
      (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) + AA*QQ))*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
      AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
      Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
      Power(AA,4)*Power(BB,2)*CC*X6))) +
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (-((Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(LL*X1) + AA*X5)) +
      (-(EE*LL) + AA*PP)*(Power(AA,5)*Power(BB,3)*CC*DD*LL*X1 -
      Power(AA,6)*Power(BB,2)*CC*DD*MM*X1 +
      Power(AA,6)*Power(BB,2)*CC*DD*MM*X2 -
      Power(AA,6)*Power(BB,3)*DD*NN*X2 +
      Power(AA,6)*Power(BB,3)*DD*NN*X3 -
      Power(AA,6)*Power(BB,3)*CC*OO*X3 +
      Power(AA,6)*Power(BB,3)*CC*OO*X4 -
      Power(AA,6)*Power(BB,3)*CC*DD*X8)))))/ ((-(EE*GG*II*LL*MM*NN) -
      BB*FF*HH*LL*NN*PP + BB*FF*II*LL*NN*PP + AA*FF*HH*MM*NN*PP -
      AA*FF*II*MM*NN*PP + AA*GG*II*MM*NN*PP + BB*EE*HH*LL*NN*QQ -
      BB*EE*II*LL*NN*QQ - AA*EE*HH*MM*NN*QQ + AA*EE*II*MM*NN*QQ +
      CC*EE*HH*LL*MM*RR - BB*EE*HH*LL*NN*RR + BB*EE*II*LL*NN*RR -
      AA*CC*HH*MM*PP*RR + AA*BB*HH*NN*PP*RR - AA*BB*II*NN*PP*RR -
      CC*EE*GG*LL*MM*SS + BB*EE*GG*LL*NN*SS + BB*CC*FF*LL*PP*SS -
      AA*CC*FF*MM*PP*SS + AA*CC*GG*MM*PP*SS - AA*BB*GG*NN*PP*SS -
      BB*CC*EE*LL*QQ*SS + AA*CC*EE*MM*QQ*SS + CC*EE*GG*LL*MM*TT -
      BB*CC*FF*LL*PP*TT + AA*CC*FF*MM*PP*TT - AA*CC*GG*MM*PP*TT +
      BB*CC*EE*LL*QQ*TT - AA*CC*EE*MM*QQ*TT - BB*CC*EE*LL*RR*TT +
      AA*BB*CC*PP*RR*TT)* (-((-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* ((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
      Power(AA,6)*Power(BB,3)*DD*HH*NN +
      Power(AA,6)*Power(BB,3)*DD*II*NN -
      Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR))) - (-(EE*LL) +
      AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
      Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
      Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
      (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS))* (Power(AA,4)*Power(BB,2)*JJ*NN
      - Power(AA,4)*Power(BB,2)*CC*UU)) + (-(EE*LL) + AA*PP)*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (-((-(EE*LL) +
      AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
      Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
      Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
      (Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS)) + (-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN -
      Power(AA,4)*Power(BB,2)*CC*TT))*
      (Power(AA,6)*Power(BB,3)*DD*JJ*NN -
      Power(AA,6)*Power(BB,3)*CC*JJ*OO +
      Power(AA,6)*Power(BB,3)*CC*KK*OO -
      Power(AA,6)*Power(BB,3)*CC*DD*VV)));

  d_bma = (GG*II*LL*NN*PP*X1 - CC*HH*LL*PP*RR*X1 + AA*HH*NN*PP*RR*X1 -
      AA*II*NN*PP*RR*X1 + CC*GG*LL*PP*SS*X1 - AA*GG*NN*PP*SS*X1 -
      CC*GG*LL*PP*TT*X1 + AA*CC*PP*RR*TT*X1 - FF*II*LL*NN*PP*X2 +
      EE*II*LL*NN*QQ*X2 - EE*II*LL*NN*RR*X2 + AA*II*NN*PP*RR*X2 -
      CC*FF*LL*PP*SS*X2 + AA*FF*NN*PP*SS*X2 + CC*EE*LL*QQ*SS*X2 -
      AA*EE*NN*QQ*SS*X2 + CC*FF*LL*PP*TT*X2 - CC*EE*LL*QQ*TT*X2 +
      CC*EE*LL*RR*TT*X2 - AA*CC*PP*RR*TT*X2 + FF*HH*LL*NN*PP*X3 -
      EE*HH*LL*NN*QQ*X3 + EE*HH*LL*NN*RR*X3 - AA*HH*NN*PP*RR*X3 -
      EE*GG*LL*NN*SS*X3 - AA*FF*NN*PP*SS*X3 + AA*GG*NN*PP*SS*X3 +
      AA*EE*NN*QQ*SS*X3 - EE*GG*II*LL*NN*X5 + CC*EE*HH*LL*RR*X5 -
      AA*EE*HH*NN*RR*X5 + AA*EE*II*NN*RR*X5 - CC*EE*GG*LL*SS*X5 +
      AA*EE*GG*NN*SS*X5 + CC*EE*GG*LL*TT*X5 - AA*CC*EE*RR*TT*X5 +
      EE*GG*II*LL*NN*X6 + CC*FF*HH*LL*PP*X6 - AA*FF*HH*NN*PP*X6 +
      AA*FF*II*NN*PP*X6 - AA*GG*II*NN*PP*X6 - CC*EE*HH*LL*QQ*X6 +
      AA*EE*HH*NN*QQ*X6 - AA*EE*II*NN*QQ*X6 - CC*EE*GG*LL*TT*X6 -
      AA*CC*FF*PP*TT*X6 + AA*CC*GG*PP*TT*X6 + AA*CC*EE*QQ*TT*X6 -
      CC*FF*HH*LL*PP*X7 + CC*EE*HH*LL*QQ*X7 - CC*EE*HH*LL*RR*X7 +
      AA*CC*HH*PP*RR*X7 + CC*EE*GG*LL*SS*X7 + AA*CC*FF*PP*SS*X7 -
      AA*CC*GG*PP*SS*X7 - AA*CC*EE*QQ*SS*X7)/ (EE*GG*II*LL*MM*NN +
      BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP - AA*FF*HH*MM*NN*PP +
      AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP - BB*EE*HH*LL*NN*QQ +
      BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ - AA*EE*II*MM*NN*QQ -
      CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR - BB*EE*II*LL*NN*RR +
      AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR + AA*BB*II*NN*PP*RR +
      CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS - BB*CC*FF*LL*PP*SS +
      AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS + AA*BB*GG*NN*PP*SS +
      BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS - CC*EE*GG*LL*MM*TT +
      BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT + AA*CC*GG*MM*PP*TT -
      BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT + BB*CC*EE*LL*RR*TT -
      AA*BB*CC*PP*RR*TT) + ((-(FF*HH*LL*PP) + EE*HH*LL*QQ -
      EE*HH*LL*RR + AA*HH*PP*RR + EE*GG*LL*SS + AA*FF*PP*SS -
      AA*GG*PP*SS - AA*EE*QQ*SS)*(JJ*NN - CC*UU)*
      (-(((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
      Power(AA,6)*Power(BB,3)*DD*HH*NN +
      Power(AA,6)*Power(BB,3)*DD*II*NN -
      Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR))) - (-(EE*LL) +
      AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
      Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
      Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
      (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS))*
      (-(((Power(AA,3)*Power(BB,2)*CC*FF*LL - Power(AA,4)*BB*CC*FF*MM
      + Power(AA,4)*BB*CC*GG*MM -
      Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
      (Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ))*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
      AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
      Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
      Power(AA,4)*Power(BB,2)*CC*X6))) +
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
      AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
      Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
      Power(AA,4)*Power(BB,2)*NN*X2 + Power(AA,4)*Power(BB,2)*NN*X3 -
      Power(AA,4)*Power(BB,2)*CC*X7)))) + (-((-(EE*LL) +
      AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
      Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
      Power(AA,4)*Power(BB,2)*GG*NN)* (-(EE*LL) + AA*PP) -
      (Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS)) + (-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN -
      Power(AA,4)*Power(BB,2)*CC*TT))*
      (-(((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
      Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
      Power(AA,6)*Power(BB,3)*DD*GG*NN)*(-(EE*LL) + AA*PP) -
      (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) + AA*QQ))*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
      AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
      Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
      Power(AA,4)*Power(BB,2)*CC*X6))) +
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (-((Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(LL*X1) + AA*X5)) +
      (-(EE*LL) + AA*PP)*(Power(AA,5)*Power(BB,3)*CC*DD*LL*X1 -
      Power(AA,6)*Power(BB,2)*CC*DD*MM*X1 +
      Power(AA,6)*Power(BB,2)*CC*DD*MM*X2 -
      Power(AA,6)*Power(BB,3)*DD*NN*X2 +
      Power(AA,6)*Power(BB,3)*DD*NN*X3 -
      Power(AA,6)*Power(BB,3)*CC*OO*X3 +
      Power(AA,6)*Power(BB,3)*CC*OO*X4 -
      Power(AA,6)*Power(BB,3)*CC*DD*X8)))))/ ((EE*GG*II*LL*MM*NN +
      BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP - AA*FF*HH*MM*NN*PP +
      AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP - BB*EE*HH*LL*NN*QQ +
      BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ - AA*EE*II*MM*NN*QQ -
      CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR - BB*EE*II*LL*NN*RR +
      AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR + AA*BB*II*NN*PP*RR +
      CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS - BB*CC*FF*LL*PP*SS +
      AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS + AA*BB*GG*NN*PP*SS +
      BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS - CC*EE*GG*LL*MM*TT +
      BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT + AA*CC*GG*MM*PP*TT -
      BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT + BB*CC*EE*LL*RR*TT -
      AA*BB*CC*PP*RR*TT)* (-((-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* ((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
      Power(AA,6)*Power(BB,3)*DD*HH*NN +
      Power(AA,6)*Power(BB,3)*DD*II*NN -
      Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR))) - (-(EE*LL) +
      AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
      Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
      Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
      (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
      Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS))* (Power(AA,4)*Power(BB,2)*JJ*NN
      - Power(AA,4)*Power(BB,2)*CC*UU)) + (-(EE*LL) + AA*PP)*
      (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (-((-(EE*LL) +
      AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
      Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
      Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
      (Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
      AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*CC*SS)) + (-(EE*LL) +
      AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
      Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
      AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
      Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
      Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
      Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN -
      Power(AA,4)*Power(BB,2)*CC*TT))*
      (Power(AA,6)*Power(BB,3)*DD*JJ*NN -
      Power(AA,6)*Power(BB,3)*CC*JJ*OO +
      Power(AA,6)*Power(BB,3)*CC*KK*OO -
      Power(AA,6)*Power(BB,3)*CC*DD*VV)));

  d_uma = (-X2 + X3)/CC + (GG*(-(BB*LL*PP*X1) + AA*MM*PP*X1 +
    EE*LL*MM*X2 - AA*MM*PP*X2 + BB*EE*LL*X5 - AA*EE*MM*X5 -
    BB*EE*LL*X6 + AA*BB*PP*X6))/ (CC*(EE*GG*LL*MM - BB*FF*LL*PP +
    AA*FF*MM*PP - AA*GG*MM*PP + BB*EE*LL*QQ - AA*EE*MM*QQ -
    BB*EE*LL*RR + AA*BB*PP*RR)) - (((-HH + II)/CC - (GG*(EE*LL -
    AA*PP)*(-(HH*MM) + BB*SS))/(CC*(EE*GG*LL*MM - BB*FF*LL*PP +
    AA*FF*MM*PP - AA*GG*MM*PP + BB*EE*LL*QQ - AA*EE*MM*QQ -
    BB*EE*LL*RR + AA*BB*PP*RR)))* (-(BB*GG*LL*NN*PP*X1) +
    AA*GG*MM*NN*PP*X1 + BB*CC*LL*PP*RR*X1 - AA*CC*MM*PP*RR*X1 +
    BB*FF*LL*NN*PP*X2 - AA*FF*MM*NN*PP*X2 - BB*EE*LL*NN*QQ*X2 +
    AA*EE*MM*NN*QQ*X2 - CC*EE*LL*MM*RR*X2 + BB*EE*LL*NN*RR*X2 +
    AA*CC*MM*PP*RR*X2 - AA*BB*NN*PP*RR*X2 + EE*GG*LL*MM*NN*X3 -
    BB*FF*LL*NN*PP*X3 + AA*FF*MM*NN*PP*X3 - AA*GG*MM*NN*PP*X3 +
    BB*EE*LL*NN*QQ*X3 - AA*EE*MM*NN*QQ*X3 - BB*EE*LL*NN*RR*X3 +
    AA*BB*NN*PP*RR*X3 + BB*EE*GG*LL*NN*X5 - AA*EE*GG*MM*NN*X5 -
    BB*CC*EE*LL*RR*X5 + AA*CC*EE*MM*RR*X5 + CC*EE*GG*LL*MM*X6 -
    BB*EE*GG*LL*NN*X6 - BB*CC*FF*LL*PP*X6 + AA*CC*FF*MM*PP*X6 -
    AA*CC*GG*MM*PP*X6 + AA*BB*GG*NN*PP*X6 + BB*CC*EE*LL*QQ*X6 -
    AA*CC*EE*MM*QQ*X6 - CC*EE*GG*LL*MM*X7 + BB*CC*FF*LL*PP*X7 -
    AA*CC*FF*MM*PP*X7 + AA*CC*GG*MM*PP*X7 - BB*CC*EE*LL*QQ*X7 +
    AA*CC*EE*MM*QQ*X7 + BB*CC*EE*LL*RR*X7 - AA*BB*CC*PP*RR*X7))/
    (EE*GG*II*LL*MM*NN + BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP -
    AA*FF*HH*MM*NN*PP + AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP -
    BB*EE*HH*LL*NN*QQ + BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ -
    AA*EE*II*MM*NN*QQ - CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR -
    BB*EE*II*LL*NN*RR + AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR +
    AA*BB*II*NN*PP*RR + CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS -
    BB*CC*FF*LL*PP*SS + AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS +
    AA*BB*GG*NN*PP*SS + BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS -
    CC*EE*GG*LL*MM*TT + BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT +
    AA*CC*GG*MM*PP*TT - BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT +
    BB*CC*EE*LL*RR*TT - AA*BB*CC*PP*RR*TT) - ((JJ/CC + ((EE*GG*LL*MM -
    BB*FF*LL*PP + AA*FF*MM*PP - AA*GG*MM*PP + BB*EE*LL*QQ -
    AA*EE*MM*QQ - BB*EE*LL*RR + AA*BB*PP*RR)* ((-HH + II)/CC -
    (GG*(EE*LL - AA*PP)*(-(HH*MM) + BB*SS))/(CC*(EE*GG*LL*MM -
    BB*FF*LL*PP + AA*FF*MM*PP - AA*GG*MM*PP + BB*EE*LL*QQ -
    AA*EE*MM*QQ - BB*EE*LL*RR + AA*BB*PP*RR)))* (-(JJ*NN) +
    CC*UU))/(EE*GG*II*LL*MM*NN + BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP
    - AA*FF*HH*MM*NN*PP + AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP -
    BB*EE*HH*LL*NN*QQ + BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ -
    AA*EE*II*MM*NN*QQ - CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR -
    BB*EE*II*LL*NN*RR + AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR +
    AA*BB*II*NN*PP*RR + CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS -
    BB*CC*FF*LL*PP*SS + AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS +
    AA*BB*GG*NN*PP*SS + BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS -
    CC*EE*GG*LL*MM*TT + BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT +
    AA*CC*GG*MM*PP*TT - BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT +
    BB*CC*EE*LL*RR*TT - AA*BB*CC*PP*RR*TT))*
    (-(((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
    Power(AA,6)*Power(BB,3)*DD*HH*NN +
    Power(AA,6)*Power(BB,3)*DD*II*NN -
    Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR))) - (-(EE*LL) +
    AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
    Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
    Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
    (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
    AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
    Power(AA,4)*Power(BB,2)*CC*SS))*
    (-(((Power(AA,3)*Power(BB,2)*CC*FF*LL - Power(AA,4)*BB*CC*FF*MM +
    Power(AA,4)*BB*CC*GG*MM - Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL)
    + AA*PP) - (Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ))*
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
    AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 - Power(AA,4)*BB*CC*MM*X1
    + Power(AA,4)*BB*CC*MM*X2 - Power(AA,4)*Power(BB,2)*CC*X6))) +
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR)))* (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
    AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 - Power(AA,4)*BB*CC*MM*X1
    + Power(AA,4)*BB*CC*MM*X2 - Power(AA,4)*Power(BB,2)*NN*X2 +
    Power(AA,4)*Power(BB,2)*NN*X3 - Power(AA,4)*Power(BB,2)*CC*X7))))
    + (-((-(EE*LL) + AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
    Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
    Power(AA,4)*Power(BB,2)*GG*NN)* (-(EE*LL) + AA*PP) -
    (Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
    AA*QQ))*(Power(AA,4)*BB*CC*HH*MM - Power(AA,4)*Power(BB,2)*CC*SS))
    + (-(EE*LL) + AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
    Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN -
    Power(AA,4)*Power(BB,2)*CC*TT))*
    (-(((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
    Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
    Power(AA,6)*Power(BB,3)*DD*GG*NN)*(-(EE*LL) + AA*PP) -
    (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) + AA*QQ))*
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
    AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 - Power(AA,4)*BB*CC*MM*X1
    + Power(AA,4)*BB*CC*MM*X2 - Power(AA,4)*Power(BB,2)*CC*X6))) +
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR)))* (-((Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(LL*X1) + AA*X5)) +
    (-(EE*LL) + AA*PP)*(Power(AA,5)*Power(BB,3)*CC*DD*LL*X1 -
    Power(AA,6)*Power(BB,2)*CC*DD*MM*X1 +
    Power(AA,6)*Power(BB,2)*CC*DD*MM*X2 -
    Power(AA,6)*Power(BB,3)*DD*NN*X2 +
    Power(AA,6)*Power(BB,3)*DD*NN*X3 -
    Power(AA,6)*Power(BB,3)*CC*OO*X3 +
    Power(AA,6)*Power(BB,3)*CC*OO*X4 -
    Power(AA,6)*Power(BB,3)*CC*DD*X8)))))/ (-((-(EE*LL) +
    AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR)))* ((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
    Power(AA,6)*Power(BB,3)*DD*HH*NN +
    Power(AA,6)*Power(BB,3)*DD*II*NN -
    Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR))) - (-(EE*LL) +
    AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
    Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
    Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
    (Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
    Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
    AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
    Power(AA,4)*Power(BB,2)*CC*SS))* (Power(AA,4)*Power(BB,2)*JJ*NN -
    Power(AA,4)*Power(BB,2)*CC*UU)) + (-(EE*LL) + AA*PP)*
    (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR)))* (-((-(EE*LL) +
    AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
    Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
    Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
    (Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
    AA*QQ))*(Power(AA,4)*BB*CC*HH*MM - Power(AA,4)*Power(BB,2)*CC*SS))
    + (-(EE*LL) + AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
    Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
    AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
    Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
    Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
    Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN -
    Power(AA,4)*Power(BB,2)*CC*TT))* (Power(AA,6)*Power(BB,3)*DD*JJ*NN
    - Power(AA,6)*Power(BB,3)*CC*JJ*OO +
    Power(AA,6)*Power(BB,3)*CC*KK*OO -
    Power(AA,6)*Power(BB,3)*CC*DD*VV));

  d_uba = (-X3 + X4)/DD + (II*(-(BB*GG*LL*NN*PP*X1) +
	AA*GG*MM*NN*PP*X1 + BB*CC*LL*PP*RR*X1 - AA*CC*MM*PP*RR*X1 +
	BB*FF*LL*NN*PP*X2 - AA*FF*MM*NN*PP*X2 - BB*EE*LL*NN*QQ*X2 +
	AA*EE*MM*NN*QQ*X2 - CC*EE*LL*MM*RR*X2 + BB*EE*LL*NN*RR*X2 +
	AA*CC*MM*PP*RR*X2 - AA*BB*NN*PP*RR*X2 + EE*GG*LL*MM*NN*X3 -
	BB*FF*LL*NN*PP*X3 + AA*FF*MM*NN*PP*X3 - AA*GG*MM*NN*PP*X3 +
	BB*EE*LL*NN*QQ*X3 - AA*EE*MM*NN*QQ*X3 - BB*EE*LL*NN*RR*X3 +
	AA*BB*NN*PP*RR*X3 + BB*EE*GG*LL*NN*X5 - AA*EE*GG*MM*NN*X5 -
	BB*CC*EE*LL*RR*X5 + AA*CC*EE*MM*RR*X5 + CC*EE*GG*LL*MM*X6 -
	BB*EE*GG*LL*NN*X6 - BB*CC*FF*LL*PP*X6 + AA*CC*FF*MM*PP*X6 -
	AA*CC*GG*MM*PP*X6 + AA*BB*GG*NN*PP*X6 + BB*CC*EE*LL*QQ*X6 -
	AA*CC*EE*MM*QQ*X6 - CC*EE*GG*LL*MM*X7 + BB*CC*FF*LL*PP*X7 -
	AA*CC*FF*MM*PP*X7 + AA*CC*GG*MM*PP*X7 - BB*CC*EE*LL*QQ*X7 +
	AA*CC*EE*MM*QQ*X7 + BB*CC*EE*LL*RR*X7 - AA*BB*CC*PP*RR*X7))/
	(DD*(EE*GG*II*LL*MM*NN + BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP
	- AA*FF*HH*MM*NN*PP + AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP -
	BB*EE*HH*LL*NN*QQ + BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ -
	AA*EE*II*MM*NN*QQ - CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR -
	BB*EE*II*LL*NN*RR + AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR +
	AA*BB*II*NN*PP*RR + CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS -
	BB*CC*FF*LL*PP*SS + AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS +
	AA*BB*GG*NN*PP*SS + BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS -
	CC*EE*GG*LL*MM*TT + BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT +
	AA*CC*GG*MM*PP*TT - BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT +
	BB*CC*EE*LL*RR*TT - AA*BB*CC*PP*RR*TT)) - (((-JJ + KK)/DD -
	(II*(EE*GG*LL*MM - BB*FF*LL*PP + AA*FF*MM*PP - AA*GG*MM*PP +
	BB*EE*LL*QQ - AA*EE*MM*QQ - BB*EE*LL*RR +
	AA*BB*PP*RR)*(-(JJ*NN) + CC*UU))/ (DD*(EE*GG*II*LL*MM*NN +
	BB*FF*HH*LL*NN*PP - BB*FF*II*LL*NN*PP - AA*FF*HH*MM*NN*PP +
	AA*FF*II*MM*NN*PP - AA*GG*II*MM*NN*PP - BB*EE*HH*LL*NN*QQ +
	BB*EE*II*LL*NN*QQ + AA*EE*HH*MM*NN*QQ - AA*EE*II*MM*NN*QQ -
	CC*EE*HH*LL*MM*RR + BB*EE*HH*LL*NN*RR - BB*EE*II*LL*NN*RR +
	AA*CC*HH*MM*PP*RR - AA*BB*HH*NN*PP*RR + AA*BB*II*NN*PP*RR +
	CC*EE*GG*LL*MM*SS - BB*EE*GG*LL*NN*SS - BB*CC*FF*LL*PP*SS +
	AA*CC*FF*MM*PP*SS - AA*CC*GG*MM*PP*SS + AA*BB*GG*NN*PP*SS +
	BB*CC*EE*LL*QQ*SS - AA*CC*EE*MM*QQ*SS - CC*EE*GG*LL*MM*TT +
	BB*CC*FF*LL*PP*TT - AA*CC*FF*MM*PP*TT + AA*CC*GG*MM*PP*TT -
	BB*CC*EE*LL*QQ*TT + AA*CC*EE*MM*QQ*TT + BB*CC*EE*LL*RR*TT -
	AA*BB*CC*PP*RR*TT)))* (-(((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM
	- Power(AA,6)*Power(BB,3)*DD*HH*NN +
	Power(AA,6)*Power(BB,3)*DD*II*NN -
	Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR))) - (-(EE*LL) +
	AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
	Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
	Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
	Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
	(Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
	Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
	AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
	Power(AA,4)*Power(BB,2)*CC*SS))*
	(-(((Power(AA,3)*Power(BB,2)*CC*FF*LL -
	Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
	Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
	(Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ))*
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
	AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
	Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
	Power(AA,4)*Power(BB,2)*CC*X6))) +
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR)))* (-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
	AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
	Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
	Power(AA,4)*Power(BB,2)*NN*X2 + Power(AA,4)*Power(BB,2)*NN*X3
	- Power(AA,4)*Power(BB,2)*CC*X7)))) + (-((-(EE*LL) +
	AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
	Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
	Power(AA,4)*Power(BB,2)*GG*NN)* (-(EE*LL) + AA*PP) -
	(Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
	AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
	Power(AA,4)*Power(BB,2)*CC*SS)) + (-(EE*LL) +
	AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
	Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN
	- Power(AA,4)*Power(BB,2)*CC*TT))*
	(-(((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
	Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
	Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
	Power(AA,6)*Power(BB,3)*DD*GG*NN)*(-(EE*LL) + AA*PP) -
	(Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
	Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) + AA*QQ))*
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(LL*X1) + AA*X5)) + (-(EE*LL) +
	AA*PP)*(Power(AA,3)*Power(BB,2)*CC*LL*X1 -
	Power(AA,4)*BB*CC*MM*X1 + Power(AA,4)*BB*CC*MM*X2 -
	Power(AA,4)*Power(BB,2)*CC*X6))) +
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR)))* (-((Power(AA,5)*Power(BB,3)*CC*DD*EE*LL
	- Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(LL*X1) + AA*X5)) +
	(-(EE*LL) + AA*PP)*(Power(AA,5)*Power(BB,3)*CC*DD*LL*X1 -
	Power(AA,6)*Power(BB,2)*CC*DD*MM*X1 +
	Power(AA,6)*Power(BB,2)*CC*DD*MM*X2 -
	Power(AA,6)*Power(BB,3)*DD*NN*X2 +
	Power(AA,6)*Power(BB,3)*DD*NN*X3 -
	Power(AA,6)*Power(BB,3)*CC*OO*X3 +
	Power(AA,6)*Power(BB,3)*CC*OO*X4 -
	Power(AA,6)*Power(BB,3)*CC*DD*X8)))))/ (-((-(EE*LL) +
	AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR)))* ((Power(AA,6)*Power(BB,2)*CC*DD*HH*MM -
	Power(AA,6)*Power(BB,3)*DD*HH*NN +
	Power(AA,6)*Power(BB,3)*DD*II*NN -
	Power(AA,6)*Power(BB,3)*CC*II*OO)*(-(EE*LL) + AA*PP)*
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR))) - (-(EE*LL) +
	AA*PP)*((Power(AA,5)*Power(BB,3)*CC*DD*FF*LL -
	Power(AA,6)*Power(BB,2)*CC*DD*FF*MM +
	Power(AA,6)*Power(BB,2)*CC*DD*GG*MM -
	Power(AA,6)*Power(BB,3)*DD*GG*NN)* (-(EE*LL) + AA*PP) -
	(Power(AA,5)*Power(BB,3)*CC*DD*EE*LL -
	Power(AA,6)*Power(BB,2)*CC*DD*EE*MM)*(-(FF*LL) +
	AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
	Power(AA,4)*Power(BB,2)*CC*SS))*
	(Power(AA,4)*Power(BB,2)*JJ*NN -
	Power(AA,4)*Power(BB,2)*CC*UU)) + (-(EE*LL) + AA*PP)*
	(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR)))* (-((-(EE*LL) +
	AA*PP)*((Power(AA,3)*Power(BB,2)*CC*FF*LL -
	Power(AA,4)*BB*CC*FF*MM + Power(AA,4)*BB*CC*GG*MM -
	Power(AA,4)*Power(BB,2)*GG*NN)*(-(EE*LL) + AA*PP) -
	(Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) +
	AA*QQ))*(Power(AA,4)*BB*CC*HH*MM -
	Power(AA,4)*Power(BB,2)*CC*SS)) + (-(EE*LL) +
	AA*PP)*(-((Power(AA,3)*Power(BB,2)*CC*EE*LL -
	Power(AA,4)*BB*CC*EE*MM)*(-(FF*LL) + AA*QQ)) + (-(EE*LL) +
	AA*PP)*(Power(AA,4)*BB*CC*GG*MM -
	Power(AA,2)*BB*CC*(-(AA*BB*FF*LL) + Power(AA,2)*FF*MM +
	Power(AA,2)*BB*RR)))* (Power(AA,4)*BB*CC*HH*MM -
	Power(AA,4)*Power(BB,2)*HH*NN + Power(AA,4)*Power(BB,2)*II*NN
	- Power(AA,4)*Power(BB,2)*CC*TT))*
	(Power(AA,6)*Power(BB,3)*DD*JJ*NN -
	Power(AA,6)*Power(BB,3)*CC*JJ*OO +
	Power(AA,6)*Power(BB,3)*CC*KK*OO -
	Power(AA,6)*Power(BB,3)*CC*DD*VV));

  // printf("specs: nby: %g, nmy: %g, ty: %g, umy: %g, uby: %g\n", get_bby(), get_bmy(), get_ty(), get_umy(), get_uby());

  // if (am_avoiding_sub_MT and (get_bmy() < 0.0 or get_ty() < 0.0 or get_umy() < 0.0 or get_uby() < 0.0)) {
   if (am_avoiding_sub_MT and (get_bmy() < 0.0 or get_ty() < 0.0 or get_umy() < 0.0)) {
    // printf("Onebound domain under MT, retrying...\n");
    // fprintf(stderr, "Onebound domain under MT, retrying...\n");
    return false;
  }
  return true;
}

double Dynein_onebound::get_binding_rate() {
  if (abs(get_uby()) < MICROTUBULE_BINDING_DISTANCE) {
    if (binding_mode == GIBBS_FULL) {
      if (am_debugging_conversions) printf("Creating bothbound from onebound to test energy\n");
      double dG_spring = Dynein_bothbound(this, rand, true).get_PE() - get_PE();
      if (isnan(dG_spring)) return 0.0;
      return low_affinity_binding_rate * exp(-dG_spring/kb/T);
    }
    // shouldn't this at least depend on the height above the microtubule
    else if (binding_mode == EXPONENTIAL_UNBINDING) {
      return low_affinity_binding_rate;
    }
    else if (binding_mode == GIBBS_BD) {
      if (am_debugging_conversions) printf("Creating bothbound from onebound to test energy\n");
      double dG_spring_BD;
      double bb_binding_equilibrium = bothbound_pre_powerstroke_internal_angles.nba;
      if (state == NEARBOUND)
	dG_spring_BD = Power((Dynein_bothbound(this, rand, true).get_fba() - bb_binding_equilibrium), 2)*cb/2.0;
      else
	dG_spring_BD = Power((Dynein_bothbound(this, rand, true).get_nba() - bb_binding_equilibrium), 2)*cb/2.0;
      if (isnan(dG_spring_BD)) return 0.0;
      return low_affinity_binding_rate * exp(-dG_spring_BD/kb/T);
    }
  }
  return 0.0;
}


//NOTE: this turns off the ability for dynein to fall off of the microtubule
double Dynein_onebound::get_unbinding_rate() {
  //return high_affinity_unbinding_preexponential_factor*exp(-PE_bba/kb/T);
  return 0;
}

/*** Set positions, velocities and forces ***/
void Dynein_onebound::set_bbx(double d) {   // onebound
  bbx = d;
}

void Dynein_onebound::set_bby(double d) {   // onebound
  bby = d;
}

void Dynein_onebound::set_bba(double d) {   // onebound
  bba = d;
}

void Dynein_onebound::set_bma(double d) {   // onebound
  bma = d;
}

void Dynein_onebound::set_uma(double d) {   // onebound
  uma = d;
}

void Dynein_onebound::set_uba(double d) {   // onebound
  uba = d;
}

/*** Get coordinates ***/

double Dynein_onebound::get_bbx() {
  return bbx;
}

double Dynein_onebound::get_bby(){
  return bby;
}

double Dynein_onebound::get_bmx() {
  return Ls * cos(bba) + bbx;
}

double Dynein_onebound::get_bmy() {
  return Ls * sin(bba) + bby;
}

double Dynein_onebound::get_tx() {
  return Lt * cos(bma) + get_bmx();
}

double Dynein_onebound::get_ty(){
  return Lt * sin(bma) + get_bmy();
}

double Dynein_onebound::get_umx() {
  return -Lt * cos(uma) + get_tx();
}

double Dynein_onebound::get_umy(){
  return -Lt * sin(uma) + get_ty();
}

double Dynein_onebound::get_ubx() {
  return -Ls * cos(uba) + get_umx();
}

double Dynein_onebound::get_uby(){
  return -Ls * sin(uba) + get_umy();
}

/*** Get Cartesian Velocities ***/

double Dynein_onebound::get_d_bbx() {
  return 0;
}

double Dynein_onebound::get_d_bmx() {
  return Ls * d_bba * -sin(bba);
}

double Dynein_onebound::get_d_tx() {
  return Lt * d_bma * -sin(bma) + get_d_bmx();
}

double Dynein_onebound::get_d_umx() {
  return Lt * d_uma * sin(uma) + get_d_tx();
}

double Dynein_onebound::get_d_ubx() {
  return Ls * d_uba * sin(uba) + get_d_umx();
}

double Dynein_onebound::get_d_bby() {
  return 0;
}

double Dynein_onebound::get_d_bmy() {
    return Ls * d_bba * cos(bba);
}

double Dynein_onebound::get_d_ty() {
    return Lt * d_bma * cos(bma) + get_d_bmy();
}

double Dynein_onebound::get_d_umy() {
    return Lt * d_uma * -cos(uma) + get_d_ty();
}

double Dynein_onebound::get_d_uby() {
  return Ls * d_uba * -cos(uba) + get_d_umy();
}
