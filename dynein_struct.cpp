#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"

/* *********************************** DYNEIN FUNCTIONS ****************************************** */

Dynein::Dynein(double bba_init, double bma_init, double fma_init, double fba_init,
               State s, forces *internal_test, forces *brownian_test, equilibrium_angles* eq_angles) {
  bbx = 0;
  bby = 0;
  mode = PRE_POWERSTROKE;

  rand = new MTRand(RAND_INIT_SEED);

  bba = bba_init;
  bma = bma_init;
  fma = fma_init;
  fba = fba_init;
  
  state = s;
  internal_testcase = internal_test;
  brownian_testcase = brownian_test;
  
  if (eq_angles) {
     eq = *eq_angles; 
  } else if (state == BOTHBOUND) {
    eq = bothbound_pre_powerstroke_internal_angles;
  } else if (state == NEARBOUND or state == FARBOUND) {
    eq = near_farbound_post_powerstroke_internal_angles;
  }
  
  update_velocities();
}

void Dynein::update_brownian_forces() {
  if (brownian_testcase) {
    r = *brownian_testcase; // just copy over forces!
  } else {
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.bbx, &r.bby);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.bmx, &r.bmy);
    rand->gauss2(sqrt(2*kb*T/(gt*dt)), &r.tx, &r.ty);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.fmx, &r.fmy);
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.fbx, &r.fby);
  }
}

void Dynein::update_internal_forces() {
  if (internal_testcase) {
    f = *internal_testcase;
  } else {
    f.bbx = 0;     f.bby = 0;   // Initialize forces to zero
    f.bmx = 0;     f.bmy = 0;
    f.tx  = 0;     f.ty  = 0;
    f.fmx = 0;     f.fmy = 0;
    f.fbx = 0;     f.fby = 0;

    double T, f1, f2, f1x, f1y, f2x, f2y;
    
    T = cb*(bba - eq.bba);
    // printf("binding angle from equilibrium: %g\n", bba - eq.bba);
    f2 = T/ls;
    f2x = f2 * sin(bba);
    f2y = f2 * -cos(bba);
    f.bmx += f2x;
    f.bmy += f2y;
    f.bbx += -f2x; // Equal and opposite forces!  :)
    f.bby += -f2y; // Equal and opposite forces!  :)

    T = cm*((bma + M_PI - bba) - eq.ba);
    // printf("bound motor from equilibrium: %g\n", (bma + M_PI - bba) - eq.ba);
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

    T = ct*(-((fma - bma) - eq.ta));
    // printf("tail from equilibrium: %g\n", -((fma - bma) - eq.ta));
    f1 = T / lt;
    f2 = T / lt;
    f1x = f1 * sin(bba);
    f1y = f1 * -cos(bba);
    f2x = f2 * sin(fma);
    f2y = f2 * -cos(fma);
    f.bmx += f1x;
    f.bmy += f1y;
    f.fmx += f2x;
    f.fmy += f2y;
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = cm*((fma + M_PI - fba) - eq.fa);
    // printf("free motor from equilibrium: %g\n", (fma + M_PI - fba) - eq.fa);
    f1 = T / lt;
    f2 = T / ls;
    f1x = f1 * sin(fma);
    f1y = f1 * -cos(fma);
    f2x = f2 * sin(fba);
    f2y = f2 * -cos(fba);
    f.tx  += f1x;
    f.ty  += f1y;
    f.fbx += f2x;
    f.fby += f2y;
    f.fmx += -(f1x + f2x);
    f.fmy += -(f1y + f2y);
    
    if (get_bmy() < 0) f.bmy += MICROTUBULE_REPULSION_FORCE * fabs(get_bmy());
    if (get_ty()  < 0) f.ty  += MICROTUBULE_REPULSION_FORCE * fabs(get_ty());
    if (get_fmy() < 0) f.fmy += MICROTUBULE_REPULSION_FORCE * fabs(get_fmy());
    if (get_fby() < 0) f.fby += MICROTUBULE_REPULSION_FORCE * fabs(get_fby());
  }
}

void Dynein::set_state(State s) {
  state = s;
}

void Dynein::update_velocities() {
  update_internal_forces();
  update_brownian_forces();

  if (state == NEARBOUND || state == FARBOUND) {
    update_velocities_onebound();
  } else if (state == BOTHBOUND) {
    update_velocities_bothbound();
  }
}

void Dynein::switch_to_bothbound() {
  // At this time, actually just switch to near/farbound. Eventually implement bothbound.
  double temp_bba = bba;
  double temp_bma = bma;
  double temp_fma = fma;
  double temp_fba = fba;

  bbx = get_fbx();
  bby = 0;
  
  bba = temp_fba;
  bma = temp_fma;
  fma = temp_bma;
  fba = temp_bba;

  if (state == NEARBOUND) state = FARBOUND;
  else if (state == FARBOUND) state = NEARBOUND;
}

void Dynein::unbind() {
  state = UNBOUND;
}

void Dynein::update_velocities_onebound() {
  float A1, A2, A3, A4;  // Start solving for velocities with matrix solution in derivation.pdf
  float B1, B2, B3, B4;
  float C1, C2, C3, C4;
  float D1, D2, D3, D4;
  float X1, X2, X3, X4;
  float Nbb, Nml, Nmr, Nbr;
  float D;
  
  A1 = -4*ls;
  A2 = -3*lt*(sin(bma)*sin(bba) + cos(bma)*cos(bba));
  A3 = 2*lt*(sin(fma)*sin(bba) + cos(fma)*cos(bba));
  A4 = ls*(sin(fba)*sin(bba) + cos(fba)*cos(bba));
  B1 = -3*ls*(sin(bba)*sin(bma) + cos(bba)*cos(bma));
  B2 = -3*lt;
  B3 = +2*lt*(sin(fma)*sin(bma) + cos(fma)*cos(bma));
  B4 = + ls*(sin(fba)*sin(bma) + cos(fba)*cos(bma));
  C1 = - 2*ls*(sin(bba)*sin(fma) + cos(bba)*cos(fma));
  C2 = - 2*lt*(sin(bma)*sin(fma) + cos(bma)*cos(fma));
  C3 = 2*lt;
  C4 = ls*(sin(fba)*sin(fma) + cos(fba)*cos(fma));
  D1 = ls*(cos(fba)*cos(bba) + sin(fba)*sin(bba));
  D2 = lt*(cos(fba)*cos(bma) + sin(fba)*sin(bma));
  D3 = -lt*(cos(fba)*cos(fma) + sin(fba)*sin(fma));
  D4 = -ls;

  X1 = (- 1/gm*f.bmy - 1/gt*f.ty - 1/gm*f.fmy - 1/gb*f.fby - r.bmy - r.ty - r.fmy - r.fby)*cos(bba)
      + ( 1/gm*f.bmx + 1/gt*f.tx + 1/gm*f.fmx + 1/gb*f.fbx + r.bmx + r.tx + r.fmx + r.fbx )*sin(bba);
  
  X2 = (-1/gt*f.ty - 1/gm*f.fmy - 1/gb*f.fby - r.ty - r.fmy - r.fby)*cos(bma) + (1/gt*f.tx + 1/gm*f.fmx + 1/gb*f.fbx + r.tx + r.fmx + r.fbx)*sin(bma);

  X3 = (-r.fmy -r.fby - 1/gm*f.fmy - 1/gb*f.fby)*cos(fma) + (r.fmx + r.fbx + 1/gm*f.fmx + 1/gb*f.fbx)*sin(fma);
  
  X4 = (r.fby + 1/gb*f.fby)*cos(fba) - (r.fbx + 1/gb*f.fbx)*sin(fba);

  Nbb = (-B2*C4*D3*X1 + B2*C3*D4*X1 + A4*C3*D2*X2 - A3*C4*D2*X2 - A4*C2*D3*X2 + A2*C4*D3*X2 + A3*C2*D4*X2 - A2*C3*D4*X2 + A4*B2*D3*X3 - A3*B2*D4*X3 - A4*B2*C3*X4 + A3*B2*C4*X4
	   + B4*(-C3*D2*X1 + C2*D3*X1 + A3*D2*X3 - A2*D3*X3 - A3*C2*X4 + A2*C3*X4) + B3*(C4*D2*X1 - C2*D4*X1 - A4*D2*X3 + A2*D4*X3 + A4*C2*X4 - A2*C4*X4));

  Nml = (B1*C4*D3*X1 - B1*C3*D4*X1 - A4*C3*D1*X2 + A3*C4*D1*X2 + A4*C1*D3*X2 - A1*C4*D3*X2 - A3*C1*D4*X2 + A1*C3*D4*X2 - A4*B1*D3*X3 + A3*B1*D4*X3 + A4*B1*C3*X4 - A3*B1*C4*X4 
	   + B4*(C3*D1*X1 - C1*D3*X1 - A3*D1*X3 + A1*D3*X3 + A3*C1*X4 - A1*C3*X4) + B3*(-C4*D1*X1 + C1*D4*X1 + A4*D1*X3 - A1*D4*X3 - A4*C1*X4 + A1*C4*X4));

  Nmr = (-B1*C4*D2*X1 + B1*C2*D4*X1 + A4*C2*D1*X2 - A2*C4*D1*X2 - A4*C1*D2*X2 + A1*C4*D2*X2 + A2*C1*D4*X2 - A1*C2*D4*X2 + A4*B1*D2*X3 - A2*B1*D4*X3 - A4*B1*C2*X4 + A2*B1*C4*X4
	   + B4*(-C2*D1*X1 + C1*D2*X1 + A2*D1*X3 - A1*D2*X3 - A2*C1*X4 + A1*C2*X4) + B2*(C4*D1*X1 - C1*D4*X1 - A4*D1*X3 + A1*D4*X3 + A4*C1*X4 - A1*C4*X4));

  Nbr = (B1*C3*D2*X1 - B1*C2*D3*X1 - A3*C2*D1*X2 + A2*C3*D1*X2 + A3*C1*D2*X2 - A1*C3*D2*X2 - A2*C1*D3*X2 + A1*C2*D3*X2 - A3*B1*D2*X3 + A2*B1*D3*X3 + A3*B1*C2*X4 - A2*B1*C3*X4 
	   + B3*(C2*D1*X1 - C1*D2*X1 - A2*D1*X3 + A1*D2*X3 + A2*C1*X4 - A1*C2*X4) + B2*(-C3*D1*X1 + C1*D3*X1 + A3*D1*X3 - A1*D3*X3 - A3*C1*X4 + A1*C3*X4));

  D = A2*B4*C3*D1 - A2*B3*C4*D1 - A1*B4*C3*D2 + A1*B3*C4*D2 - A2*B4*C1*D3 + A1*B4*C2*D3 + A2*B1*C4*D3 - A1*B2*C4*D3 + A4*(B3*C2*D1 - B2*C3*D1 - B3*C1*D2 + B1*C3*D2 + B2*C1*D3 
	   - B1*C2*D3)+ A2*B3*C1*D4 - A1*B3*C2*D4 - A2*B1*C3*D4 + A1*B2*C3*D4 + A3*(-B4*C2*D1 + B2*C4*D1 + B4*C1*D2 - B1*C4*D2 - B2*C1*D4 + B1*C2*D4);
  
  assert(D != 0);

  d_bba = Nbb/D;
  d_bma = Nml/D;
  d_fma = Nmr/D;
  d_fba = Nbr/D;
}

void Dynein::update_velocities_bothbound() {
  // To be implemented
}

double Dynein::get_binding_probability() {
  if (get_fby() < MICROTUBULE_BINDING_DISTANCE) {
    return 0.01;
  } else {
    return 0;
  }
}

double Dynein::get_unbinding_probability() {
  if (f.bby + r.bby >= UNBINDING_FORCE) { // bad, doesn't take into account forces on other domains
    printf("unbinding...\n");
    return 1.0;
  } else return 0.0;
}


/*** Set positions, velocities and forces ***/
void Dynein::set_bbx(double d) {
  bbx = d;
}

void Dynein::set_bby(double d) {
  bby = d;
}

void Dynein::set_bba(double d) {
  bba = d;
}

void Dynein::set_bma(double d) {
  bma = d;
}

void Dynein::set_fma(double d) {
  fma = d;
}

void Dynein::set_fba(double d) {
  fba = d;
}

/*** Angular Velocities ***/

double Dynein::get_d_bba() {
  return d_bba;
}

double Dynein::get_d_bma() {
  return d_bma;
}

double Dynein::get_d_fma() {
  return d_fma;
}

double Dynein::get_d_fba() {
  return d_fba;
}

/*** Get coordinates ***/

double Dynein::get_bbx() {
  return bbx;
}

double Dynein::get_bby(){
  return bby;
}

double Dynein::get_bmx() {
  return ls * cos(get_bba()) + bbx;
}

double Dynein::get_bmy(){
  return ls * sin(get_bba()) + bby;
}

double Dynein::get_tx() {
  return ls * cos(get_bba()) + lt * cos(get_bma()) + bbx;
}

double Dynein::get_ty(){
  return ls * sin(get_bba()) + lt * sin(get_bma()) + bby;
}

double Dynein::get_fmx() {
  return ls * cos(get_bba()) + lt * cos(get_bma()) - lt * cos(get_fma()) + bbx;
}

double Dynein::get_fmy(){
  return ls * sin(get_bba()) + lt * sin(get_bma()) - lt * sin(get_fma()) + bby;
}

double Dynein::get_fbx() {
  return ls * cos(get_bba()) + lt * cos(get_bma()) - lt * cos(get_fma()) - ls * cos(get_fba()) + bbx;
}

double Dynein::get_fby(){
  return ls * sin(get_bba()) + lt * sin(get_bma()) - lt * sin(get_fma()) - ls * sin(get_fba()) + bby;
}

/*** Get Cartesian Velocities ***/

double Dynein::get_d_bbx() {
  return 0;
}

double Dynein::get_d_bmx() {
  return ls * d_bba * -sin(bba);
}

double Dynein::get_d_tx() {
  return lt * d_bma * -sin(bma) + get_d_bmx();
}

double Dynein::get_d_fmx() {
  return lt * d_fma * sin(fma) + get_d_tx();
}

double Dynein::get_d_fbx() {
  return ls * d_fba * sin(fba) + get_d_fmx();
}

double Dynein::get_d_bby() {
  return 0;
}

double Dynein::get_d_bmy() {
  return ls * d_bba * cos(bba);
}

double Dynein::get_d_ty() {
  return lt * d_bma * cos(bma) + get_d_bmy();
}

double Dynein::get_d_fmy() {
  return lt * d_fma * -cos(fma) + get_d_ty();
}

double Dynein::get_d_fby() {
  return ls * d_fba * -cos(fba) + get_d_fmy();
}

/*** Get forces ***/
forces Dynein::get_internal() {
  return f;
}

forces Dynein::get_brownian() {
  return r;
}

/*** Get angles ***/

double Dynein::get_bba() {
  return bba;
}

double Dynein::get_bma() {
  return bma;
}

double Dynein::get_fma() {
  return fma;
}

double Dynein::get_fba() {
  return fba;
}

State Dynein::get_state() {
  return state;
}

/*** Get energies ***/

double Dynein::get_PE() {
  return 0;
}

double Dynein::get_KE() {
  return 0;
}

void Dynein::log(double t, FILE* data_file) {
  fprintf(data_file, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",
          get_KE(), get_PE(), get_KE() + get_PE(), t, get_bbx(), get_bby(), get_bmx(), get_bmy(),
          get_tx(), get_ty(), get_fmx(), get_fmy(), get_fbx(), get_fby(), state);
}

void Dynein::resetLog() {
	FILE* data_file = fopen("data.txt", "w");
	FILE* config_file = fopen("config.txt", "w");
	
	fprintf(config_file, "#gb\tgm\tgt\tdt\truntime?\tstate\n");
	fprintf(config_file, "%g\t%g\t%g\t%g\t%g\t%d\n",
          (double) gb, (double) gm, (double) gt, dt, runtime, (int) state);
	fprintf(data_file,
		"#KE\tPE\t\tEnergy\t\tt\t\tbbX\tbbY\tbmx\tbmy\ttX\ttY\tfmx\tfmy\tfbx\tfby\tS\n");
	
	fclose(data_file);
	fclose(config_file);
}
