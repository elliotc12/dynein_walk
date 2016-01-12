#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"

/* ********************* ONEBOUND DYNEIN FUNCTIONS ************************** */

Dynein_onebound::Dynein_onebound(double bba_init, double bma_init,
				 double uma_init, double uba_init,
				 double bbx_init, double bby_init,
				 State s, onebound_forces *internal_test,
				 onebound_forces *brownian_test,
				 onebound_equilibrium_angles* eq_angles,
				 MTRand *r) {
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
  rand = r;
  if (!rand) rand = new MTRand();

  update_velocities();
}

void Dynein_onebound::update_brownian_forces() {
  if (brownian_testcase) {
    r = *brownian_testcase; // just copy over forces!
  } else {
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.bbx, &r.bby);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.bmx, &r.bmy);
    rand->gauss2(sqrt(2*kb*T/(gt*dt)), &r.tx, &r.ty);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.umx, &r.umy);
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.ubx, &r.uby);
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
    f2 = T/ls;
    f2x = f2 * sin(bba);
    f2y = f2 * -cos(bba);
    f.bmx += f2x;
    f.bmy += f2y;
    f.bbx += -f2x; // Equal and opposite forces!  :)
    f.bby += -f2y; // Equal and opposite forces!  :)

    T = cm*((bma + M_PI - bba) - eq.bba);
    f1 = T/ls;
    f2 = T/lt;
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

    T = ct*(uma - bma - eq.ta);  //-- this used to be the negation, this is right?
    f1 = T / lt;
    f2 = T / lt;
    f1x = f1 * sin(bma);
    f1y = f1 * -cos(bma);
    f2x = f2 * -sin(uma);  // not sure if these angles are right?
    f2y = f2 * cos(uma);
    f.bmx += f1x;
    f.bmy += f1y;
    f.umx += f2x;
    f.umy += f2y;
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = cm*((uma + M_PI - uba) - eq.uma);
    f1 = T / lt;
    f2 = T / ls;
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

    if (get_bmy() < 0) f.bmy += MICROTUBULE_REPULSION_FORCE * fabs(get_bmy());
    if (get_ty()  < 0) f.ty  += MICROTUBULE_REPULSION_FORCE * fabs(get_ty());
    if (get_umy() < 0) f.umy += MICROTUBULE_REPULSION_FORCE * fabs(get_umy());
    if (get_uby() < 0) f.uby += MICROTUBULE_REPULSION_FORCE * fabs(get_uby());
  }
}

// void Dynein_onebound::switch_to_bothbound() {
//   double temp_nma;
//   double temp_fma;

//   if (state == NEARBOUND) {
//     temp_nma = M_PI + bma - bba;
//     temp_fma = M_PI + uma - uba;
//     nbx = get_bbx();
//     L = fbx - nbx;
//   } else {
//     temp_nma = M_PI + uma - uba;
//     temp_fma = M_PI + bma - bba;
//     nbx = get_fbx();
//     L = bbx - fbx;
//   }

//   distance_traveled += fabs(get_ubx() - get_bbx());
//   steps++;

//   state = BOTHBOUND;
//   nma = temp_nma;
//   fma = temp_fma;
// }

// void Dynein_onebound::switch_to_nearbound() {
//   //nearbound -> bma is nma
//   double temp_bba = get_nba();
//   double temp_bma = nma;
//   double temp_uma = fma;
//   double temp_uba = get_fba();

//   bbx = get_nbx();
//   bby = 0;

//   bba = temp_bba;
//   bma = temp_bma;
//   uma = temp_uma;
//   uba = temp_uba;

//   state = NEARBOUND;
// }

// void Dynein_onebound::switch_to_farbound() {
//   //nearbound -> bma is fma
//   double temp_bba = get_fba();
//   double temp_bma = fma;
//   double temp_uma = nma;
//   double temp_uba = get_nba();

//   bbx = get_fbx();
//   bby = 0;

//   bba = temp_bba;
//   bma = temp_bma;
//   uma = temp_uma;
//   uba = temp_uba;

//   state = FARBOUND;
// }

void Dynein_onebound::update_velocities() {
  update_internal_forces();
  update_brownian_forces();

  float A1, A2, A3, A4;  // Start solving for velocities with matrix solution in derivation.pdf
  float B1, B2, B3, B4;
  float C1, C2, C3, C4;
  float D1, D2, D3, D4;
  float X1, X2, X3, X4;
  float Nbb, Nml, Nmr, Nbr;
  float D;

  A1 = -4*ls;
  A2 = -3*lt*(sin(bma)*sin(bba) + cos(bma)*cos(bba));
  A3 = 2*lt*(sin(uma)*sin(bba) + cos(uma)*cos(bba));
  A4 = ls*(sin(uba)*sin(bba) + cos(uba)*cos(bba));
  B1 = -3*ls*(sin(bba)*sin(bma) + cos(bba)*cos(bma));
  B2 = -3*lt;
  B3 = +2*lt*(sin(uma)*sin(bma) + cos(uma)*cos(bma));
  B4 = +ls*(sin(uba)*sin(bma) + cos(uba)*cos(bma));
  C1 = -2*ls*(sin(bba)*sin(uma) + cos(bba)*cos(uma));
  C2 = -2*lt*(sin(bma)*sin(uma) + cos(bma)*cos(uma));
  C3 = 2*lt;
  C4 = ls*(sin(uba)*sin(uma) + cos(uba)*cos(uma));
  D1 = ls*(cos(uba)*cos(bba) + sin(uba)*sin(bba));
  D2 = lt*(cos(uba)*cos(bma) + sin(uba)*sin(bma));
  D3 = -lt*(cos(uba)*cos(uma) + sin(uba)*sin(uma));
  D4 = -ls;

  X1 = (- 1/gm*f.bmy - 1/gt*f.ty - 1/gm*f.umy - 1/gb*f.uby - r.bmy - r.ty - r.umy - r.uby)*cos(bba)
      +( 1/gm*f.bmx + 1/gt*f.tx + 1/gm*f.umx + 1/gb*f.ubx + r.bmx + r.tx + r.umx + r.ubx )*sin(bba);

  X2 = (-1/gt*f.ty - 1/gm*f.umy - 1/gb*f.uby - r.ty - r.umy - r.uby)*cos(bma)
    + (1/gt*f.tx + 1/gm*f.umx + 1/gb*f.ubx + r.tx + r.umx + r.ubx)*sin(bma);

  X3 = (-r.umy -r.uby - 1/gm*f.umy - 1/gb*f.uby)*cos(uma)
    + (r.umx + r.ubx + 1/gm*f.umx + 1/gb*f.ubx)*sin(uma);

  X4 = (r.uby + 1/gb*f.uby)*cos(uba) - (r.ubx + 1/gb*f.ubx)*sin(uba);

  Nbb = (-B2*C4*D3*X1 + B2*C3*D4*X1 + A4*C3*D2*X2 - A3*C4*D2*X2 - A4*C2*D3*X2 + A2*C4*D3*X2
	 +A3*C2*D4*X2 - A2*C3*D4*X2 + A4*B2*D3*X3 - A3*B2*D4*X3 - A4*B2*C3*X4 + A3*B2*C4*X4
	 +B4*(-C3*D2*X1 + C2*D3*X1 + A3*D2*X3 - A2*D3*X3 - A3*C2*X4 + A2*C3*X4)
	 +B3*(C4*D2*X1 - C2*D4*X1 - A4*D2*X3 + A2*D4*X3 + A4*C2*X4 - A2*C4*X4));

  Nml = (B1*C4*D3*X1 - B1*C3*D4*X1 - A4*C3*D1*X2 + A3*C4*D1*X2 + A4*C1*D3*X2 - A1*C4*D3*X2
	 -A3*C1*D4*X2 + A1*C3*D4*X2 - A4*B1*D3*X3 + A3*B1*D4*X3 + A4*B1*C3*X4 - A3*B1*C4*X4
	 +B4*(C3*D1*X1 - C1*D3*X1 - A3*D1*X3 + A1*D3*X3 + A3*C1*X4 - A1*C3*X4)
	 +B3*(-C4*D1*X1 + C1*D4*X1 + A4*D1*X3 - A1*D4*X3 - A4*C1*X4 + A1*C4*X4));

  Nmr = (-B1*C4*D2*X1 + B1*C2*D4*X1 + A4*C2*D1*X2 - A2*C4*D1*X2 - A4*C1*D2*X2 + A1*C4*D2*X2
	 +A2*C1*D4*X2 - A1*C2*D4*X2 + A4*B1*D2*X3 - A2*B1*D4*X3 - A4*B1*C2*X4 + A2*B1*C4*X4
	 +B4*(-C2*D1*X1 + C1*D2*X1 + A2*D1*X3 - A1*D2*X3 - A2*C1*X4 + A1*C2*X4)
	 +B2*(C4*D1*X1 - C1*D4*X1 - A4*D1*X3 + A1*D4*X3 + A4*C1*X4 - A1*C4*X4));

  Nbr = (B1*C3*D2*X1 - B1*C2*D3*X1 - A3*C2*D1*X2 + A2*C3*D1*X2 + A3*C1*D2*X2 - A1*C3*D2*X2
	 -A2*C1*D3*X2 + A1*C2*D3*X2 - A3*B1*D2*X3 + A2*B1*D3*X3 + A3*B1*C2*X4 - A2*B1*C3*X4
	 +B3*(C2*D1*X1 - C1*D2*X1 - A2*D1*X3 + A1*D2*X3 + A2*C1*X4 - A1*C2*X4)
	 +B2*(-C3*D1*X1 + C1*D3*X1 + A3*D1*X3 - A1*D3*X3 - A3*C1*X4 + A1*C3*X4));

  D = A2*B4*C3*D1 - A2*B3*C4*D1 - A1*B4*C3*D2 + A1*B3*C4*D2 - A2*B4*C1*D3 + A1*B4*C2*D3
         +A2*B1*C4*D3 - A1*B2*C4*D3 + A4*(B3*C2*D1 - B2*C3*D1 - B3*C1*D2 + B1*C3*D2 + B2*C1*D3
         -B1*C2*D3)+ A2*B3*C1*D4 - A1*B3*C2*D4 - A2*B1*C3*D4 + A1*B2*C3*D4
         +A3*(-B4*C2*D1 + B2*C4*D1 + B4*C1*D2 - B1*C4*D2 - B2*C1*D4 + B1*C2*D4);

  assert(D != 0);

  d_bba = Nbb/D;
  d_bma = Nml/D;
  d_uma = Nmr/D;
  d_uba = Nbr/D;
}

double Dynein_onebound::get_binding_rate() {
  if (get_uby() < MICROTUBULE_BINDING_DISTANCE) {
    double bound_energy = 0.5*cb*square(fba - eq.bba);
    return 1e10*exp(-bound_energy/kb/T); // per second
  } else {
    return 0;
  }
}

double Dynein_onebound::get_unbinding_rate() {
  if (f.bby + r.bby >= UNBINDING_FORCE) { // bad, doesn't take into account forces on other domains?
    printf("unbinding...\n");
    return 1.0;
  } else return 0.0;
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

void Dynein_onebound::set_nma(double d) {   // bothbound
  assert(C == BOTHBOUND);
  nma = d;
}

/*** Angular Velocities ***/

double Dynein_onebound::get_d_bba() {
  return d_bba;
}

double Dynein_onebound::get_d_bma() {
  return d_bma;
}

double Dynein_onebound::get_d_fma() {
  return d_fma;
}

double Dynein_onebound::get_d_fba() {
  return d_fba;
}

double Dynein_onebound::get_d_bma() {  // bothbound
  int pm = (bma > M_PI) ? -1 : 1; // sign of d_nma depends on value of nma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - ln*ln) / (2*lt*ls),2))
    * (ln / (lt*ls)) * d_ln;
}

double Dynein_onebound::get_d_uma() {  // bothbound
  int pm = (uma > M_PI) ? -1 : 1; // sign of d_fma depends on value of fma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - lf*lf) / (2*lt*ls),2))
    * (lf / (lt*ls)) * d_lf;
}

/*** Get coordinates ***/

double Dynein_onebound::get_bbx() {
  return bbx;
}

double Dynein_onebound::get_bby(){
  return bby;
}

double Dynein_onebound::get_bmx() {
  return ls * cos(get_bba()) + bbx;
}

double Dynein_onebound::get_bmy(){
  return ls * sin(get_bba()) + bby;
}

double Dynein_onebound::get_tx() {
  return ls * cos(get_bba()) + lt * cos(get_bma()) + bbx;
}

double Dynein_onebound::get_ty(){
  return ls * sin(get_bba()) + lt * sin(get_bma()) + bby;
}

double Dynein_onebound::get_fmx() {
  return ls * cos(get_bba()) + lt * cos(get_bma()) - lt * cos(get_fma()) + bbx;
}

double Dynein_onebound::get_umy(){
  return ls * sin(get_bba()) + lt * sin(get_bma()) - lt * sin(get_fma()) + bby;
}

double Dynein_onebound::get_fbx() {
  return ls * cos(get_bba()) + lt * cos(get_bma())
    - lt * cos(get_fma()) - ls * cos(get_fba()) + bbx;
}

double Dynein_onebound::get_uby(){
  return ls * sin(get_bba()) + lt * sin(get_bma())
    - lt * sin(get_fma()) - ls * sin(get_fba()) + bby;
}

double Dynein_onebound::get_nbx() {   // bothbound
  return nbx;
}

double Dynein_onebound::get_nmx() {   // bothbound
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*(L*ln))) p
         + pm*acos((pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) / (2*(ln*ls))));
}

double Dynein_onebound::get_tx() {   // bothbound
  return (pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*L);
}

double Dynein_onebound::get_fmx() {   // bothbound
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(lf, 2) - pow(ln, 2)) / (2*(L*lf)))
         + pm*acos((pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) / (2*(lf*ls))));;
}

double Dynein_onebound::get_fbx() {   // bothbound
  return L;
}

double Dynein_onebound::get_nby() {   // bothbound
  return 0;
}

double Dynein_onebound::get_nmy() {   // bothbound
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*(L*ln)))
         + pm*acos((pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) / (2*(ln*ls))));
}

double Dynein_onebound::get_ty() {   // bothbound
  return ln*sqrt(1 - pow(pow(L, 2) + pow(ln, 2) - pow(lf, 2), 2)
  / (4*(pow(L, 2)*pow(ln, 2))));
}

double Dynein_onebound::get_fmy() {   // bothbound
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(lf, 2) - pow(ln, 2)) / (2*(L*lf)))
         + pm*acos((pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) / (2*(lf*ls))));
}

double Dynein_onebound::get_fby() {   // bothbound
  return 0;
}


/*** Get Cartesian Velocities ***/

double Dynein_onebound::get_d_bbx() {
  return 0;
}

double Dynein_onebound::get_d_bmx() {
  return ls * d_bba * -sin(bba);
}

double Dynein_onebound::get_d_tx() {
  return lt * d_bma * -sin(bma) + get_d_bmx();
}

double Dynein_onebound::get_d_fmx() {
  return lt * d_fma * sin(fma) + get_d_tx();
}

double Dynein_onebound::get_d_fbx() {
  return ls * d_fba * sin(fba) + get_d_fmx();
}

double Dynein_onebound::get_d_bby() {
  return 0;
}

double Dynein_onebound::get_d_bmy() {
    return ls * d_bba * cos(bba);
}

double Dynein_onebound::get_d_ty() {
    return lt * d_bma * cos(bma) + get_d_bmy();
}

double Dynein_onebound::get_d_umy() {
    return lt * d_fma * -cos(fma) + get_d_ty();
}

double Dynein_onebound::get_d_uby() {
  return ls * d_fba * -cos(fba) + get_d_umy();
}

/*** Get forces ***/
onebound_forces Dynein_onebound::get_internal() {
  return f;
}

onebound_forces Dynein_onebound::get_brownian() {
  return r;
}

/*** Get angles ***/

double Dynein_onebound::get_bba() {
  return bba;
}

double Dynein_onebound::get_bma() {
  return bma;
}

double Dynein_onebound::get_fma() {
  return fma;
}

double Dynein_onebound::get_fba() {
  return fba;
}

// State Dynein_onebound::get_state() {
//   return state;
// }

/*** Get energies ***/

double Dynein_onebound::get_PE() {
  return 0;
}

double Dynein_onebound::get_KE() {
  return 0;
}

void Dynein_onebound::log(double t, FILE* data_file) {
  fprintf(data_file, "%.2g\t%.2g\t%.2g\t%.5g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t"
	  "%.4f\t%.4f\t%d\n",
          get_KE(), get_PE(), get_KE() + get_PE(), t, get_bbx(), get_bby(), get_bmx(), get_bmy(),
          get_tx(), get_ty(), get_fmx(), get_umy(), get_fbx(), get_uby(), state);
}

void Dynein_onebound::log_run(float runtime) {
  FILE* data_file = fopen("run_data.txt", "w");

  float run_length = get_bbx();
  float ave_step_dist = distance_traveled / steps;
  float ave_step_time = runtime / steps;

  printf("\n\n***********Run data**********\n");
  printf("Run length: %f nm\n", run_length);
  printf("Distance traveled: %f nm\n", distance_traveled);
  printf("Steps: %d\n", steps);
  printf("Average step length: %f nm\n", ave_step_dist);
  printf("Average step time: %g s\n\n\n", ave_step_time);
  fprintf(data_file, "Run length \tDistance traveled \tSteps \tAve step length \tAve step time\n");
  fprintf(data_file, "%f\t%f\t%d\t%f\t%g\n",
	  run_length, distance_traveled, steps, ave_step_dist, ave_step_time);
  fclose(data_file);
}

void Dynein_onebound::resetLog() {
  FILE* data_file = fopen("data.txt", "w");
  FILE* config_file = fopen("config.txt", "w");

  fprintf(config_file, "#gb\tgm\tgt\tdt\truntime?\tstate\n");
  fprintf(config_file, "%g\t%g\t%g\t%g\t%g\t%d\n",
          (double) gb, (double) gm, (double) gt, dt, runtime, (int) state);
  fprintf(data_file,
	  "#KE\tPE\tEnergy\tt\tnbX\tnbY\tnmx\tnmy\ttX\ttY\tfmx\tfmy\tfbx\tuby\tS\n");

  fclose(data_file);
  fclose(config_file);
}
