#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

double bla_t;
double mla_t;
double mra_t;
double bra_t;

double d_bla_t;
double d_mla_t;
double d_mra_t;
double d_bra_t;

double blx_t = 0;
double bly_t = 0;

double r_blx_t = 0;    double  r_bly_t = 0;
double r_mlx_t = 0;    double  r_mly_t = 0;
double r_tx_t  = 0;	   double  r_ty_t  = 0;
double r_mrx_t = 0;    double  r_mry_t = 0;
double r_brx_t = 0;    double  r_bry_t = 0;

double f_blx_t = 0;    double  f_bly_t = 0;
double f_mlx_t = 0;    double  f_mly_t = 0;
double f_tx_t  = 1;    double  f_ty_t  = 1;
double f_mrx_t = 1;    double  f_mry_t = 1;
double f_brx_t = 1;    double  f_bry_t = 1;

const int state_t = LEFTBOUND;

/*** Prototypes ***/
double get_d_blx();
double get_d_mlx();
double get_d_tx(); 
double get_d_mrx();
double get_d_brx();

double get_d_bly();
double get_d_mly();
double get_d_ty(); 
double get_d_mry();
double get_d_bry();

double get_blx();
double get_mlx();
double get_mrx();
double get_brx();

double get_bly();
double get_mly();
double get_mry();
double get_bry();


/*** Get coordinate velocities ***/

double get_d_blx() {
	if (state_t == LEFTBOUND) return 0;
	else if (state_t == BOTHBOUND) return 0;
	else return ls * d_bla_t * sin(-bla_t) + get_d_mlx();
}

double get_d_mlx() {
	if (state_t == LEFTBOUND) return ls * d_bla_t * -sin(bla_t);
	else if (state_t == BOTHBOUND) return -d_bla_t * ls * sin(bla_t);
	else return lt * d_mla_t * sin(-mla_t) + get_d_tx();
}

double get_d_tx() {
	if (state_t == LEFTBOUND) return lt * d_mla_t * -sin(mla_t) + get_d_mlx();
	else if (state_t == BOTHBOUND) return lt/2 * (-d_mra_t * sin(mra_t) + -d_mla_t * sin(mla_t) + get_d_mlx() + get_d_mrx());
	else return lt * d_mra_t * -sin(mra_t) + get_d_mrx();
}

double get_d_mrx() {
	if (state_t == LEFTBOUND) return lt * d_mra_t * sin(-mra_t) + get_d_tx();
	else if (state_t == BOTHBOUND) return -d_bra_t * ls * sin(bra_t);
	else return ls * d_bra_t * -sin(bra_t);
}

double get_d_brx() {
	if (state_t == LEFTBOUND) return ls * d_bra_t * sin(-bra_t) + get_d_mrx();
	else if (state_t == BOTHBOUND) return 0;
	else return 0;
}

double get_d_bly() {
	if (state_t == LEFTBOUND) return 0;
	else if (state_t == BOTHBOUND) return 0;
	else return ls * -d_bla_t * cos(-bla_t) + get_d_mly();
}

double get_d_mly() {
	if (state_t == LEFTBOUND) return ls * d_bla_t * cos(bla_t);
	else if (state_t == BOTHBOUND) return d_bla_t * ls * cos(bla_t);
	else return lt * -d_mla_t * cos(-mla_t) + get_d_ty();
}

double get_d_ty() {
	if (state_t == LEFTBOUND) return lt * d_mla_t * cos(mla_t) + get_d_mly();
	else if (state_t == BOTHBOUND) return lt/2 * (d_mra_t * cos(mra_t) + d_mla_t * cos(mla_t) + get_d_mly() + get_d_mry());
	else return lt * d_mra_t * cos(mra_t) + get_d_mry();
}

double get_d_mry() {
	if (state_t == LEFTBOUND) return lt * -d_mra_t * cos(-mra_t) + get_d_ty();
	else if (state_t == BOTHBOUND) return d_bra_t * ls * cos(mra_t);
	else return ls * d_bra_t * cos(bra_t);
}

double get_d_bry() {
	if (state_t == LEFTBOUND) return ls * -d_bra_t * cos(-bra_t) + get_d_mry();
	else if (state_t == BOTHBOUND) return 0;
	else return 0;
}

/*** Get coordinates ***/

double get_blx() {
	if (state_t == LEFTBOUND) return blx_t;
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_bly(){
	if (state_t == LEFTBOUND) return bly_t;
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_mlx() {
	if (state_t == LEFTBOUND) return ls * cos(bla_t) + get_blx();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_mly(){
	if (state_t == LEFTBOUND) return ls * sin(bla_t) + get_bly();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_tx() {
	if (state_t == LEFTBOUND) return ls * cos(bla_t) + lt * cos(mla_t) + get_blx();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_ty(){
	if (state_t == LEFTBOUND) return ls * sin(bla_t) + lt * sin(mla_t) + get_bly();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_mrx() {
	if (state_t == LEFTBOUND) return ls * cos(bla_t) + lt * cos(mla_t) + lt * cos(-mra_t) + get_blx();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_mry(){
	if (state_t == LEFTBOUND) return ls * sin(bla_t) + lt * sin(mla_t) + lt * sin(-mra_t) + get_bly();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_brx() {
	if (state_t == LEFTBOUND) return ls * cos(bla_t) + lt * cos(mla_t) + lt * cos(-mra_t) + ls * cos(-bra_t) + get_blx();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}

double get_bry(){
	if (state_t == LEFTBOUND) return ls * sin(bla_t) + lt * sin(mla_t) + lt * sin(-mra_t) + ls * sin(-bra_t) + get_bly();
	else if (state_t == RIGHTBOUND) return 0;
	else return 0;
}


int main() {
	
	bla_t = (108.0 / 180) * M_PI;
	mla_t = (36.0 / 180) * M_PI;
	mra_t = (36.0 / 180) * M_PI;
	bra_t = (108.0 / 180) * M_PI;
	
	Dynein* dyn = new Dynein(bla_t, mla_t, mra_t, bra_t);

	d_bla_t = dyn->get_d_bla();
	d_mla_t = dyn->get_d_mla();
	d_mra_t = dyn->get_d_mra();
	d_bra_t = dyn->get_d_bra();
	
	printf("one: %f\n", (get_d_mrx() + lt*sin(bra_t-M_PI)*d_bra_t - (1/g)*f_brx_t - r_brx_t ) * 1/(get_brx() - get_mrx()));
	printf("two: %f\n", (get_d_mry() - lt*cos(bra_t-M_PI)*d_bra_t - (1/g)*f_bry_t - r_bry_t ) * 1/(get_bry() - get_mry()));
	
	exit(EXIT_SUCCESS);
}













