#include <stdlib.h>
#include <fstream>

#include "dynein_struct.h"


/* *********************************** DYNEIN FUNCTIONS ****************************************** */

Dynein::Dynein(double bla_init, double mla_init, double mra_init, double bra_init) {
	
	read_init_file();
	
	blx = 0;
	bly = 0;
	
	d_bla = 0;
	d_mla = 0;
	d_mra = 0;
	d_bra = 0;
	
	bla = bla_init;
	mla = mla_init;
	mra = mra_init;
	bra = bra_init;
	
	state = LEFTBOUND;
	
}

void Dynein::read_init_file() {
	// Eventually put initialization parameter reading in here
}

void Dynein::set_state(states s) {
	state = s;
}

void Dynein::update_protein() {
	
	r_blx = 0;     r_bly = 0;
	r_mlx = 0;     r_mly = 0;
	r_mrx = 0;     r_mry = 0;
	r_brx = 0;     r_bry = 0;

	f_blx = 0;     f_bly = 0;
	f_mlx = 0;     f_mly = 0;
	f_mrx = 0;     f_mry = 0;
	f_brx = 0;     f_bry = 0;
	
	int A1, A2, A3, A4;
	int B1, B2, B3, B4;
	int C1, C2, C3, C4;
	int D1, D2, D3, D4;
	int X1, X2, X3, X4;
	int Nbl, Nml, Nmr, Nbr;
	int D;
	
	if (state == LEFTBOUND) {

			A1 = -3*sin(bla)*sin(bla) - 3*cos(bla)*cos(bla) - 1;
			A2 = (-3*ls*sin(mla)*sin(bla) + -3*ls*cos(mla)*cos(bla))/lt;
			A3 = (ls*sin(mra-M_PI)*sin(bra-M_PI) - ls*cos(mra-M_PI)*cos(bra-M_PI))/lt;
			A4 = -1;
			B1 = (2*lt*cos(bla)*cos(mra-M_PI) + -2*lt*sin(bla)*sin(mra-M_PI))/ls;
			B2 = (-2*sin(mla)*sin(mra-M_PI) + -2*cos(mra-M_PI)*cos(mla));
			B3 = (-cos(mra-M_PI)*cos(mra-M_PI) + -sin(mra-M_PI)*sin(mra-M_PI) - 1);
			B4 = (-lt*sin(bra-M_PI)*sin(mra-M_PI) + -lt*cos(bra-M_PI)*cos(mra-M_PI))/ls;
			C1 = (-3*lt*sin(bla)*sin(mla) + -3*lt*cos(bla)*cos(mla))/ls;
			C2 = (-2*sin(mla)*sin(mla) + -2*cos(mla)*cos(mla) - 1);
			C3 = (-2*sin(mra-M_PI)*sin(mla) + -2*cos(mra-M_PI)*cos(mla));
			C4 = (-lt*sin(bra-M_PI)*sin(mla) + -lt*cos(bra-M_PI)*cos(mla))/ls;
			D1 = (-3*sin(bla)*sin(bla) - 3*cos(bla)*cos(bla) -1);
			D2 = (-3*ls*sin(mla)*sin(bla) + -3*ls*cos(mla)*cos(bla))/lt;
			D3 = (-2*ls*sin(mra-M_PI)*sin(bla) + -2*ls*cos(mra-M_PI)*cos(bla))/lt;
			D4 = (-sin(bla)*sin(bra-M_PI) + -cos(bra-M_PI)*cos(bla));


			X_1 = (sin(bra-M_PI)*f_brx + cos(bra-M_PI)*f_bry)/(g*lt) + (sin(bra-M_PI)*r_brx + cos(bra-M_PI)*r_bry)/(lt);
			X_2 = -(sin(mra-M_PI)f_mrx + sin(mra-M_PI)f_brx - cos(mra-M_PI)*f_mry - cos(mra-M_PI)*f_bry)/(g*ls) + 
		-(sin(mra-M_PI)*r_mrx + sin(mra-M_PI)*r_brx - cos(mra-M_PI)*r_mry - cos(mra-M_PI)*r_bry)/(ls);
			X_3 = -(sin(mla)*f_tx + -cos(mla)*f_yt + sin(mla)f_mrx + -cos(mla)*f_mry + sin(mla)f_brx + -cos(mla)*f_bry)(g*ls) + 
		-(sin(mla)*r_tx + -cos(mla)r_ty + sin(mla)*r_mrx + -cos(mla)*r_mry + sin(mla)*r_brx + -cos(mla)*r_bry)/(ls);
			X_4 = -(sin(bla)*f_mlx + -cos(bla)*f_mly + sin(bla)F_{xt } + -cos(bla)*f_ty + sin(bla)*f_mrx + -cos(bla)*f_mry + sin(bla)*f_brx + -cos(bla)*f_bry)/(g*lt)
		 + -(sin(bla)*r_mlx +-cos(bla)*r_mly + sin(bla)*r_tx + -cos(bla)*r_ty + sin(bla)*r_mrx + -cos(bla)*r_mry + sin(bla)*r_brx + -cos(bla)*r_bry)/(lt);

			Nbl =
		(-B_2*C_4*D_3*X_1 + B_2*C_3*D_4*X_1 + A_4*C_3*D_2*X_2 - A_3*C_4*D_2*X_2 - A_4*C_2*D_3*X_2 + A_2*C_4*D_3*X_2 +
		A_3*C_2*D_4*X_2 - A_2*C_3*D_4*X_2 + A_4*B_2*D_3*X_3 - A_3*B_2*D_4*X_3 - A_4*B_2*C_3*X_4 + A_3*B_2*C_4*X_4 +
		B_4*(-C_3*D_2*X_1 + C_2*D_3*X_1 + A_3*D_2*X_3 - A_2*D_3*X_3 - A_3*C_2*X_4 + A_2*C_3*X_4) + B_3*(C_4*D_2*X_1 - 
		C_2*D_4*X_1 - A_4*D_2*X_3 + A_2*D_4*X_3 + A_4*C_2*X_4 - A_2*C_4*X_4));


			Nml =
		(B_1*C_4*D_3*X_1 - B_1*C_3*D_4*X_1 - A_4*C_3*D1*X_2 + A_3*C_4*D1*X_2 + A_4*C_1*D_3*X_2 - A_1*C_4*D_3*X_2 - 
		A_3*C_1*D_4*X_2 + A_1*C_3*D_4*X_2 - A_4*B_1*D_3*X_3 + A_3*B_1*D_4*X_3 + A_4*B_1*C_3*X_4 - A_3*B_1*C_4*X_4 +
		B_4*(C_3*D1*X_1 - C_1*D_3*X_1 - A_3*D1*X_3 + A_1*D_3*X_3 + A_3*C_1*X_4 - A_1*C_3*X_4) + B_3*(-C_4*D1*X_1 + 
		C_1*D_4*X_1 + A_4*D1*X_3 - A_1*D_4*X_3 - A_4*C_1*X_4 + A_1*C_4*X_4));


			Nmr =
		(-B_1*C_4*D_2*X_1 + B_1*C_2*D_4*X_1 + A_4*C_2*D1*X_2 - A_2*C_4*D1*X_2 - A_4*C_1*D_2*X_2 + A_1*C_4*D_2*X_2 + 
		A_2*C_1*D_4*X_2 - A_1*C_2*D_4*X_2 + A_4*B_1*D_2*X_3 - A_2*B_1*D_4*X_3 - A_4*B_1*C_2*X_4 + A_2*B_1*C_4*X_4 +
		B_4*(-C_2*D1*X_1 + C_1*D_2*X_1 + A_2*D1*X_3 - A_1*D_2*X_3 - A_2*C_1*X_4 + A_1*C_2*X_4) + 
		B_2*(C_4*D1*X_1 - C_1*D_4*X_1 - A_4*D1*X_3 + A_1*D_4*X_3 + A_4*C_1*X_4 - A_1*C_4*X_4));


			Nbr =
		(B_1*C_3*D_2*X_1 - B_1*C_2*D_3*X_1 - A_3*C_2*D1*X_2 + A_2*C_3*D1*X_2 + A_3*C_1*D_2*X_2 - A_1*C_3*D_2*X_2 - 
		A_2*C_1*D_3*X_2 + A_1*C_2*D_3*X_2 - A_3*B_1*D_2*X_3 + A_2*B_1*D_3*X_3 + A_3*B_1*C_2*X_4 - A_2*B_1*C_3*X_4 +
		B_3*(C_2*D1*X_1 - C_1*D_2*X_1 - A_2*D1*X_3 + A_1*D_2*X_3 + A_2*C_1*X_4 - A_1*C_2*X_4) + 
		B_2*(-C_3*D1*X_1 + C_1*D_3*X_1 + A_3*D1*X_3 - A_1*D_3*X_3 - A_3*C_1*X_4 + A_1*C_3*X_4));

			D = A_2*B_4*C_3*D1 - A_2*B_3*C_4*D1 - A_1*B_4*C_3*D_2 + A_1*B_3*C_4*D_2 - A_2*B_4*C_1*D_3 + A_1*B_4*C_2*D_3 + 
		A_2*B_1*C_4*D_3 - A_1*B_2*C_4*D_3 + A_4*(B_3*C_2*D1 - B_2*C_3*D1 - B_3*C_1*D_2 + B_1*C_3*D_2 + B_2*C_1*D_3 - B_1*C_2*D_3) + 
		A_2*B_3*C_1*D_4 - A_1*B_3*C_2*D_4 - A_2*B_1*C_3*D_4 + A_1*B_2*C_3*D_4 + A_3*(-B_4*C_2*D1 + B_2*C_4*D1 + B_4*C_1*D_2 - B_1*C_4*D_2 - B_2*C_1*D_4 + B_1*C_2*D_4);
	
		d_bla = Nbl/D;
		d_mla = Nml/D;
		d_mra = Nmr/D;
		d_bra = Nbr/D;
	} 
}

/*** Set positions and velocities ***/

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

void Dynein::set_d_bla(double d) {
	d_bla = d;
}

void Dynein::set_d_mla(double d) {
	d_mla = d;
}

void Dynein::set_d_mra(double d) {
	d_mra = d;
}

void Dynein::set_d_bra(double d) {
	d_bra = d;
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
	if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + lt * cos(-get_mra()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else return 0;
}

double Dynein::get_mry(){
	if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + lt * sin(-get_mra()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else return 0;
}

double Dynein::get_brx() {
	if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + lt * cos(-get_mra()) + ls * cos(-get_bra()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else return 0;
}

double Dynein::get_bry(){
	if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + lt * sin(-get_mra()) + ls * sin(-get_bra()) + bly;
	else if (state == RIGHTBOUND) return 0;
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

states Dynein::get_state() {
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
	//fprintf(data_file, "%.6f\t%12.6f\t%12.6f\t%.3f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%d\n", 
	//get_KE(), get_PE(), get_KE() + get_PE(), t, blx, bly, mlx, mly, tx, ty, mrx, mry, brx, bry, get_state());
	fclose(data_file);
}
