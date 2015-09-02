#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"
extern int runtime;

/* *********************************** DYNEIN FUNCTIONS ****************************************** */

Dynein::Dynein(double bla_init, double mla_init, double mra_init, double bra_init,
               State s, forces *internal_test, forces *brownian_test, equilibrium_angles* eq_angles) {
  blx = 0;
  bly = 0;
  mode = PRE_POWERSTROKE;
  
  bla = bla_init;
  mla = mla_init;
  mra = mra_init;
  bra = bra_init;
  
  state = s;
  internal_testcase = internal_test;
  brownian_testcase = brownian_test;
  
  if (eq_angles) {
     eq = *eq_angles; 
  } else {
    eq = pre_powerstroke_internal_angles;
  }
  
  update_velocities();
}

void Dynein::update_brownian_forces() {
  if (brownian_testcase) {
    r = *brownian_testcase; // just copy over forces!
  } else {
    rand.gauss2(sqrt(2*kb*T/gb), &r.blx, &r.bly);
    rand.gauss2(sqrt(2*kb*T/gm), &r.mlx, &r.mly);
    rand.gauss2(sqrt(2*kb*T/gt), &r.tx, &r.ty);
    rand.gauss2(sqrt(2*kb*T/gm), &r.mrx, &r.mry);
    rand.gauss2(sqrt(2*kb*T/gb), &r.brx, &r.bry);
  }
}

void Dynein::update_internal_forces() {
  if (internal_testcase) {
    f = *internal_testcase;
  } else {
    f.blx = 0;     f.bly = 0;
    f.mlx = 0;     f.mly = 0;
    f.tx  = 0;     f.ty  = 0;
    f.mrx = 0;     f.mry = 0;
    f.brx = 0;     f.bry = 0;

    double T, f1, f2, f1x, f1y, f2x, f2y;
    
    T = bla - eq.bla;
    f2 = T/ls;
    f2x = f2 * sin(bla);
    f2y = f2 * -cos(bla);
    f.mlx += f2x;
    f.mly += f2y;
    f.blx += -f2x; // Equal and opposite forces!  :)
    f.bly += -f2y; // Equal and opposite forces!  :)

    T = (mla + M_PI - bla) - eq.la;
    f1 = T/ls;
    f2 = T/lt;
    f1x = f1 * sin(bla);
    f1y = f1 * -cos(bla);
    f2x = f2 * sin(mla);
    f2y = f2 * -cos(mla);
    f.blx += f1x;
    f.bly += f1y;
    f.tx  += f2x;
    f.ty  += f2y;
    f.mlx += -(f1x + f2x);
    f.mly += -(f1y + f2y);

    T = -((mra - mla) - eq.ta);
    f1 = T / lt;
    f2 = T / lt;
    f1x = f1 * sin(bla);
    f1y = f1 * -cos(bla);
    f2x = f2 * sin(mra);
    f2y = f2 * -cos(mra);
    f.mlx += f1x;
    f.mly += f1y;
    f.mrx += f2x;
    f.mry += f2y;
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = ((mra + M_PI - bra) - eq.ra);
    f1 = T / lt;
    f2 = T / ls;
    f1x = f1 * sin(mra);
    f1y = f1 * -cos(mra);
    f2x = f2 * sin(bra);
    f2y = f2 * -cos(bra);
    f.tx  += f1x;
    f.ty  += f1y;
    f.brx += f2x;
    f.bry += f2y;
    f.mrx += -(f1x + f2x);
    f.mry += -(f1y + f2y);
  }
}

void Dynein::set_state(State s) {
  state = s;
}

void Dynein::update_velocities() {
  
  float A1, A2, A3, A4;
  float B1, B2, B3, B4;
  float C1, C2, C3, C4;
  float D1, D2, D3, D4;
  float X1, X2, X3, X4;
  float Nbl, Nml, Nmr, Nbr;
  float D;

  update_internal_forces();
  update_brownian_forces();
  
    if (state == LEFTBOUND) {

    A1 = -4*ls;
    A2 = -3*lt*(sin(mla)*sin(bla) + cos(mla)*cos(bla));
    A3 = 2*lt*(sin(mra)*sin(bla) + cos(mra)*cos(bla));
    A4 = ls*(sin(bra)*sin(bla) + cos(bra)*cos(bla));
    B1 = -3*ls*(sin(bla)*sin(mla) + cos(bla)*cos(mla));
    B2 = -3*lt;
    B3 = +2*lt*(sin(mra)*sin(mla) + cos(mra)*cos(mla));
    B4 = + ls*(sin(bra)*sin(mla) + cos(bra)*cos(mla));
    C1 = - 2*ls*(sin(bla)*sin(mra) + cos(bla)*cos(mra));
    C2 = - 2*lt*(sin(mla)*sin(mra) + cos(mla)*cos(mra));
    C3 = 2*lt;
    C4 = ls*(sin(bra)*sin(mra) + cos(bra)*cos(mra));
    D1 = ls*(cos(bra)*cos(bla) + sin(bra)*sin(bla));
    D2 = lt*(cos(bra)*cos(mla) + sin(bra)*sin(mla));
    D3 = -lt*(cos(bra)*cos(mra) + sin(bra)*sin(mra));
    D4 = -ls;

    X1 = (- 1/gm*f.mly - 1/gt*f.ty - 1/gm*f.mry - 1/gb*f.bry - r.mly - r.ty - r.mry - r.bry)*cos(bla)
      + ( 1/gm*f.mlx + 1/gt*f.tx + 1/gm*f.mrx + 1/gb*f.brx + r.mlx + r.tx + r.mrx + r.brx )*sin(bla);
  
    X2 = (-1/gt*f.ty - 1/gm*f.mry - 1/gb*f.bry - r.ty - r.mry - r.bry)*cos(mla) + (1/gt*f.tx + 1/gm*f.mrx + 1/gb*f.brx + r.tx + r.mrx + r.brx)*sin(mla);

    X3 = (-r.mry -r.bry - 1/gm*f.mry - 1/gb*f.bry)*cos(mra) + (r.mrx + r.brx + 1/gm*f.mrx + 1/gb*f.brx)*sin(mra);
  
    X4 = (r.bry + 1/gb*f.bry)*cos(bra) - (r.brx + 1/gb*f.brx)*sin(bra);

    Nbl = (-B2*C4*D3*X1 + B2*C3*D4*X1 + A4*C3*D2*X2 - A3*C4*D2*X2 - A4*C2*D3*X2 + A2*C4*D3*X2 + A3*C2*D4*X2 - A2*C3*D4*X2 + A4*B2*D3*X3 - A3*B2*D4*X3 - A4*B2*C3*X4 + A3*B2*C4*X4
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

    d_bla = Nbl/D;
    d_mla = Nml/D;
    d_mra = Nmr/D;
    d_bra = Nbr/D;
  } 
}

/*** Set positions, velocities and forces ***/

void Dynein::set_blx(double d) {
  blx = d;
}

void Dynein::set_bly(double d) {
  bly = d;
}

void Dynein::set_bla(double d) {
  bla = d;
}

void Dynein::set_mla(double d) {
  mla = d;
}

void Dynein::set_mra(double d) {
  mra = d;
}

void Dynein::set_bra(double d) {
  bra = d;
}

/*** Angular Velocities ***/

double Dynein::get_d_bla() {
  return d_bla;
}

double Dynein::get_d_mla() {
  return d_mla;
}

double Dynein::get_d_mra() {
  return d_mra;
}

double Dynein::get_d_bra() {
  return d_bra;
}

/*** Get coordinates ***/

double Dynein::get_blx() {
  if (state == LEFTBOUND) return blx;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_bly(){
  if (state == LEFTBOUND) return bly;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_mlx() {
  if (state == LEFTBOUND) return ls * cos(get_bla()) + blx;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_mly(){
  if (state == LEFTBOUND) return ls * sin(get_bla()) + bly;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_tx() {
  if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + blx;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_ty(){
  if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + bly;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_mrx() {
  if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + -lt * cos(get_mra()) + blx;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_mry(){
  if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + -lt * sin(get_mra()) + bly;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_brx() {
  if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + -lt * cos(get_mra()) + -ls * cos(get_bra()) + blx;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

double Dynein::get_bry(){
  if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + -lt * sin(get_mra()) + -ls * sin(get_bra()) + bly;
  else if (state == RIGHTBOUND) return 0;
  else return 0;
}

/*** Get Cartesian Velocities ***/

double Dynein::get_d_blx() {
  if (state == LEFTBOUND) return 0;
  else if (state == BOTHBOUND) return 0;
  else return ls * d_bla * sin(-bla) + get_d_mlx();
}

double Dynein::get_d_mlx() {
  if (state == LEFTBOUND) return ls * d_bla * -sin(bla);
  else if (state == BOTHBOUND) return -d_bla * ls * sin(bla);
  else return lt * d_mla * sin(-mla) + get_d_tx();
}

double Dynein::get_d_tx() {
  if (state == LEFTBOUND) return lt * d_mla * -sin(mla) + get_d_mlx();
  else if (state == BOTHBOUND) return lt/2 * (-d_mra * sin(mra) + -d_mla * sin(mla) + get_d_mlx() + get_d_mrx());
  else return lt * d_mra * -sin(mra) + get_d_mrx();
}

double Dynein::get_d_mrx() {
  if (state == LEFTBOUND) return lt * d_mra * sin(mra) + get_d_tx();
  else if (state == BOTHBOUND) return -d_bra * ls * sin(bra);
  else return ls * d_bra * -sin(bra);
}

double Dynein::get_d_brx() {
  if (state == LEFTBOUND) return ls * d_bra * sin(bra) + get_d_mrx();
  else if (state == BOTHBOUND) return 0;
  else return 0;
}

double Dynein::get_d_bly() {
  if (state == LEFTBOUND) return 0;
  else if (state == BOTHBOUND) return 0;
  else return ls * -d_bla * cos(-bla) + get_d_mly();
}

double Dynein::get_d_mly() {
  if (state == LEFTBOUND) return ls * d_bla * cos(bla);
  else if (state == BOTHBOUND) return d_bla * ls * cos(bla);
  else return lt * -d_mla * cos(-mla) + get_d_ty();
}

double Dynein::get_d_ty() {
  if (state == LEFTBOUND) return lt * d_mla * cos(mla) + get_d_mly();
  else if (state == BOTHBOUND) return lt/2 * (d_mra * cos(mra) + d_mla * cos(mla) + get_d_mly() + get_d_mry());
  else return lt * d_mra * cos(mra) + get_d_mry();
}

double Dynein::get_d_mry() {
  if (state == LEFTBOUND) return lt * d_mra * -cos(mra) + get_d_ty();
  else if (state == BOTHBOUND) return d_bra * ls * cos(mra);
  else return ls * d_bra * cos(bra);
}

double Dynein::get_d_bry() {
  if (state == LEFTBOUND) return ls * d_bra * -cos(bra) + get_d_mry();
  else if (state == BOTHBOUND) return 0;
  else return 0;
}


/*** Get forces ***/
forces Dynein::get_internal() {
  return f;
}

/*** Get Brownian forces ***/

forces Dynein::get_brownian() {
  return r;
}

/*** Get angles ***/

double Dynein::get_bla() {
  return bla;
}

double Dynein::get_mla() {
  return mla;
}

double Dynein::get_mra() {
  return mra;
}

double Dynein::get_bra() {
  return bra;
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

void Dynein::log(double t) {
  FILE* data_file = fopen("data.txt", "a+");
  fprintf(data_file, "%.6f\t%12.6f\t%12.6f\t%.3f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%d\n", 
  get_KE(), get_PE(), get_KE() + get_PE(), t, get_blx(), get_bly(), get_mlx(), get_mly(), get_tx(), get_ty(), get_mrx(), get_mry(), get_brx(), get_bry(), get_state());
  fclose(data_file);
}

void Dynein::resetLog() {
	FILE* data_file = fopen("data.txt", "w");
	FILE* config_file = fopen("config.txt", "w");
	
	fprintf(config_file, "#inctime\truntime\tstate\n%+.3f\t%+.3f\t%d\n", inctime, (double) runtime, (int) state);
	fprintf(data_file,
		"#KE\t\t\t\tPE\t\t\t\tEnergy\t\tt\t\tblX\t\t\tblY\t\t\tmlX\t\t\tmlY\t\t\ttX\t\t\ttY\t\t\tmrX\t\t\tmrY\t\t\tbrX\t\t\tbrY\t\t\tS\n");
	
	fclose(data_file);
	fclose(config_file);
}
