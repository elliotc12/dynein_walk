#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"

/* *********************************** DYNEIN FUNCTIONS ****************************************** */

Dynein::Dynein(double bla_init, double mla_init, double mra_init, double bra_init, State s, Mode m, Brownian_mode bm) {
  
  read_init_file();

  blx = 0;
  bly = 0;
  
  bla = bla_init;
  mla = mla_init;
  mra = mra_init;
  bra = bra_init;
  
  state = s;
  mode = m;
  bmode = bm;
  
  update_velocities();
}

void Dynein::update_brownian_forces() {
  if (bmode == TEST_NO_BROWNIAN_FORCES) {
    r_blx = 0;     r_bly = 0;
    r_mlx = 0;     r_mly = 0;
    r_tx  = 0;     r_ty  = 0;
    r_mrx = 0;     r_mry = 0;
    r_brx = 0;     r_bry = 0;
  }
  if (bmode == TEST_RIGHT_BROWNIAN_FORCES) {
    r_blx = 0;       r_bly = 0;
    r_mlx = 1.0;     r_mly = 0;
    r_tx  = 1.0;     r_ty  = 0;
    r_mrx = 1.0;     r_mry = 0;
    r_brx = 1.0;     r_bry = 0;
  }
  if (bmode == TEST_LEFT_BROWNIAN_FORCES) {
    r_blx = 0;       r_bly = 0;
    r_mlx = -1.0;    r_mly = 0;
    r_tx  = -1.0;    r_ty  = 0;
    r_mrx = -1.0;    r_mry = 0;
    r_brx = -1.0;    r_bry = 0;
  }
}

void Dynein::update_internal_forces() {
  if (mode == TEST_NO_INTERNAL_FORCES) {
    f_blx = 0;     f_bly = 0;
    f_mlx = 0;     f_mly = 0;
    f_tx  = 0;     f_ty  = 0;
    f_mrx = 0;     f_mry = 0;
    f_brx = 0;     f_bry = 0;
  }

 if (mode == TEST_LEFT_INTERNAL_FORCES) {
    f_blx = 0;       f_bly = 0;
    f_mlx = -1.0;    f_mly = 0;
    f_tx  = -1.0;    f_ty  = 0;
    f_mrx = -1.0;    f_mry = 0;
    f_brx = -1.0;    f_bry = 0;
  }
  
  if (mode == PRE_POWERSTROKE) {
    f_blx = 0;
    f_mlx = (bla - bla_0) * sin(bla);
    f_tx  = ((mla + M_PI - bla) - la_0) * sin(mla);
    f_mrx = -((mra - mla) - ta_0) * sin(mra);
    f_brx = ((mra + M_PI - bra) - ra_0) * sin(bra);

    f_bly = 0;
    f_mly = (bla - bla_0) * -cos(bla);
    f_ty  = ((mla + M_PI - bla) - la_0) * -cos(mla);
    f_mry = -((mra - mla) - ta_0) * -cos(mra);
    f_bry = ((mra + M_PI - bra) - ra_0) * -cos(bra);
  }
}

void Dynein::read_init_file() {
  // Eventually put initialization parameter reading in here
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

    X1 = (- 1/g*f_mly - 1/g*f_ty - 1/g*f_mry - 1/g*f_bry - r_mly - r_ty - r_mry - r_bry)*cos(bla)
      + ( 1/g*f_mlx + 1/g*f_tx + 1/g*f_mrx + 1/g*f_brx + r_mlx + r_tx + r_mrx + r_brx )*sin(bla);
  
    X2 = (-1/g*f_ty - 1/g*f_mry - 1/g*f_bry - r_ty - r_mry - r_bry)*cos(mla) + (1/g*f_tx + 1/g*f_mrx + 1/g*f_brx + r_tx + r_mrx + r_brx)*sin(mla);

    X3 = (-r_mry -r_bry - 1/g*f_mry - 1/g*f_bry)*cos(mra) + (r_mrx + r_brx + 1/g*f_mrx + 1/g*f_brx)*sin(mra);
  
    X4 = (r_bry + 1/g*f_bry)*cos(bra) - (r_brx + 1/g*f_brx)*sin(bra);

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
double Dynein::get_f_blx() {
  return 0;
}

double Dynein::get_f_mlx() {
  return 0;
}

double Dynein::get_f_tx() {
  return 0;
}

double Dynein::get_f_mrx() {
  return 0;
}

double Dynein::get_f_brx() {
  return 0;
}

double Dynein::get_f_bly() {
  return 0;
}

double Dynein::get_f_mly() {
  return 0;
}

double Dynein::get_f_ty() {
  return 0;
}

double Dynein::get_f_mry() {
  return 0;
}

double Dynein::get_f_bry() {
  return 0;
}

/*** Get Brownian forces ***/

double Dynein::get_r_blx() {
  return 0;
}

double Dynein::get_r_mlx() {
  return 0;
}

double Dynein::get_r_tx() {
  return 0;
}

double Dynein::get_r_mrx() {
  return 0;
}

double Dynein::get_r_brx() {
  return 0;
}

double Dynein::get_r_bly() {
  return 0;
}

double Dynein::get_r_mly() {
  return 0;
}

double Dynein::get_r_ty() {
  return 0;
}

double Dynein::get_r_mry() {
  return 0;
}

double Dynein::get_r_bry() {
  return 0;
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
