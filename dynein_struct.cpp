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

		A1 = \mn 3\sin(\Theta_{bl})\sin(\Theta_{bl})\mn + 3\cos(\Theta_{bl})\cos(\Theta_{bl}) -1 \\
		A2 = \frac{\mn 3L_{s}\sin(\Theta_{ml})\sin(\Theta_{bl}) + \mn 3L_{s}\cos(\Theta_{ml})\cos(\Theta_{bl})}{L_t} \\
		A3 = \frac{L_{s}\sin(\Theta_{mr}-\pi)\sin(\Theta_{br}-\pi) - L_{s}\cos(\Theta_{mr}-\pi)\cos(\Theta_{br}-\pi)}{L_{t}} \\
		A4 = -1 \\
		B1 = \frac{2L_{t}\cos(\Theta_{bl})\cos(\Theta_{mr}-\pi) + \mn 2L_{t}\sin(\Theta_{bl})\sin(\Theta_{mr}-\pi)}{L_{s}} \\
		B2 = \big(\mn 2\sin(\Theta_{ml})\sin(\Theta_{mr}-\pi) + \mn 2\cos(\Theta_{mr}-\pi)\cos(\Theta_{ml}) \big) \\
		B3 = \big( \mn \cos(\Theta_{mr}-\pi)\cos(\Theta_{mr}-\pi) + \mn\sin(\Theta_{mr}-\pi)\sin(\Theta_{mr}-\pi) - 1\big) \\
		B4 = \frac{\mn L_{t}\sin(\Theta_{br}-\pi)\sin(\Theta_{mr}-\pi) + \mn L_{t}\cos(\Theta_{br}-\pi)\cos(\Theta_{mr}-\pi)}{L_{s}} \\
		C1 = \frac{\mn 3L_{t}\sin(\Theta_{bl})\sin(\Theta_{ml}) + \mn3L_{t}\cos(\Theta_{bl})\cos(\Theta_{ml})}{L_s} \\
		C2 = \big(\mn 2\sin(\Theta_{ml})\sin(\Theta_{ml}) + \mn 2\cos(\Theta_{ml})\cos(\Theta_{ml}) - 1\big) \\
		C3 = \big(\mn 2\sin(\Theta_{mr}-\pi)\sin(\Theta_{ml}) + \mn 2\cos(\Theta_{mr}-\pi)\cos(\Theta_{ml}) \big) \\
		C4 = \frac{\mn L_{t}\sin(\Theta_{br}-\pi)\sin(\Theta_{ml}) + \mn L_{t}\cos(\Theta_{br}-\pi)\cos(\Theta_{ml})}{L_s} \\
		D1 = \big(\mn 3\sin(\Theta_{bl})\sin(\Theta_{bl})\mn + 3\cos(\Theta_{bl})\cos(\Theta_{bl}) -1\big) \\
		D2 = \frac{\mn 3L_{s}\sin(\Theta_{ml})\sin(\Theta_{bl}) + \mn 3L_{s}\cos(\Theta_{ml})\cos(\Theta_{bl})}{L_t} \\
		D3 = \frac{\mn 2L_{s}\sin(\Theta_{mr}-\pi)\sin(\Theta_{bl}) + \mn 2L_{s}\cos(\Theta_{mr}-\pi)\cos(\Theta_{bl})}{L_t}\dot{\Theta}_{mr} \\
		D4 = \big(\mn \sin(\Theta_{bl}) \sin(\Theta_{br}-\pi) + \mn \cos(\Theta_{br}-\pi)\cos(\Theta_{bl})\big) \\
		
		
		X_1 = \frac{\sin(\Theta_{br}-\pi)F_{xbr} + \cos(\Theta_{br}-\pi)F_{ybr}}{\gamma L_{t}} + \frac{\sin(\Theta_{br}-\pi)R_{xbr} + \cos(\Theta_{br}-\pi)R_{ybr}}{L_{t}} \\
		X_2 = \mn\frac{\sin(\Theta_{mr}-\pi)F_{xmr} + \sin(\Theta_{mr}-\pi)F_{xbr} - \cos(\Theta_{mr}-\pi)F_{ymr} - \cos(\Theta_{mr}-\pi)F_{ybr}}{L_{s}\gamma} + \\
  \mn\frac{\sin(\Theta_{mr}-\pi)R_{xmr} + \sin(\Theta_{mr}-\pi)R_{xbr} - \cos(\Theta_{mr}-\pi)R_{ymr} - \cos(\Theta_{mr}-\pi)R_{ybr}}{L_{s}}
		X_3 = \mn \frac{\sin(\Theta_{ml})F_{xt } + \mn \cos(\Theta_{ml})F_{yt } + \sin(\Theta_{ml})F_{xmr} + \mn \cos(\Theta_{ml})F_{ymr} + \sin(\Theta_{ml})F_{xbr} + \mn \cos(\Theta_{ml})F_{ybr}}{\gamma L_s} + \\
	\mn \frac{\sin(\Theta_{ml})R_{xt } + \mn \cos(\Theta_{ml})R_{yt } + \sin(\Theta_{ml})R_{xmr} + \mn \cos(\Theta_{ml})R_{ymr} + \sin(\Theta_{ml})R_{xbr} + \mn \cos(\Theta_{ml})R_{ybr}}{L_s}
		X_4 = \mn \frac{\sin(\Theta_{bl})F_{xml} + \mn\cos(\Theta_{bl})F_{yml} + \sin(\Theta_{bl})F_{xt } + \mn \cos(\Theta_{bl})F_{yt } + \sin(\Theta_{bl})F_{xmr} + \mn \cos(\Theta_{bl})F_{ymr} + \sin(\Theta_{bl})F_{xbr} + \mn \cos(\Theta_{bl})F_{ybr}}
  {\gamma L_t} \\+ \mn \frac{\sin(\Theta_{bl})R_{xml} +\mn \cos(\Theta_{bl})R_{yml} + \sin(\Theta_{bl})R_{xt } + \mn \cos(\Theta_{bl})R_{yt} + \sin(\Theta_{bl})R_{xmr} + \mn \cos(\Theta_{bl})R_{ymr} + \sin(\Theta_{bl})R_{xbr} + \mn \cos(\Theta_{bl})R_{ybr}}{L_t} \\
		
		N_{bl} =
(-B_2 C_4 D_3 X_1 + B_2 C_3 D_4 X_1 + A_4 C_3 D_2 X_2 - A_3 C_4 D_2 X_2 - A_4 C_2 D_3 X_2 + A_2 C_4 D_3 X_2 +\\
 A_3 C_2 D_4 X_2 - A_2 C_3 D_4 X_2 + A_4 B_2 D_3 X_3 - A_3 B_2 D_4 X_3 - A_4 B_2 C_3 X_4 + A_3 B_2 C_4 X_4 +\\
B_4 (-C_3 D_2 X_1 + C_2 D_3 X_1 + A_3 D_2 X_3 - A_2 D_3 X_3 - A_3 C_2 X_4 + A_2 C_3 X_4) + B_3 (C_4 D_2 X_1 - \\
C_2 D_4 X_1 - A_4 D_2 X_3 + A_2 D_4 X_3 + A_4 C_2 X_4 - A_2 C_4 X_4)) \\


		N_{ml} =
(B_1 C_4 D_3 X_1 - B_1 C_3 D_4 X_1 - A_4 C_3 D1 X_2 + A_3 C_4 D1 X_2 + A_4 C_1 D_3 X_2 - A_1 C_4 D_3 X_2 - \\
A_3 C_1 D_4 X_2 + A_1 C_3 D_4 X_2 - A_4 B_1 D_3 X_3 + A_3 B_1 D_4 X_3 + A_4 B_1 C_3 X_4 - A_3 B_1 C_4 X_4 +\\
B_4 (C_3 D1 X_1 - C_1 D_3 X_1 - A_3 D1 X_3 + A_1 D_3 X_3 + A_3 C_1 X_4 - A_1 C_3 X_4) + B_3 (-C_4 D1 X_1 + \\
C_1 D_4 X_1 + A_4 D1 X_3 - A_1 D_4 X_3 - A_4 C_1 X_4 + A_1 C_4 X_4)) \\


		N_{mr} =
(-B_1 C_4 D_2 X_1 + B_1 C_2 D_4 X_1 + A_4 C_2 D1 X_2 - A_2 C_4 D1 X_2 - A_4 C_1 D_2 X_2 + A_1 C_4 D_2 X_2 + \\
A_2 C_1 D_4 X_2 - A_1 C_2 D_4 X_2 + A_4 B_1 D_2 X_3 - A_2 B_1 D_4 X_3 - A_4 B_1 C_2 X_4 + A_2 B_1 C_4 X_4 +\\
B_4 (-C_2 D1 X_1 + C_1 D_2 X_1 + A_2 D1 X_3 - A_1 D_2 X_3 - A_2 C_1 X_4 + A_1 C_2 X_4) + \\
B_2 (C_4 D1 X_1 - C_1 D_4 X_1 - A_4 D1 X_3 + A_1 D_4 X_3 + A_4 C_1 X_4 - A_1 C_4 X_4)) \\


		N_{br} =
(B_1 C_3 D_2 X_1 - B_1 C_2 D_3 X_1 - A_3 C_2 D1 X_2 + A_2 C_3 D1 X_2 + A_3 C_1 D_2 X_2 - A_1 C_3 D_2 X_2 - \\
A_2 C_1 D_3 X_2 + A_1 C_2 D_3 X_2 - A_3 B_1 D_2 X_3 + A_2 B_1 D_3 X_3 + A_3 B_1 C_2 X_4 - A_2 B_1 C_3 X_4 +\\
B_3 (C_2 D1 X_1 - C_1 D_2 X_1 - A_2 D1 X_3 + A_1 D_2 X_3 + A_2 C_1 X_4 - A_1 C_2 X_4) + \\
B_2 (-C_3 D1 X_1 + C_1 D_3 X_1 + A_3 D1 X_3 - A_1 D_3 X_3 - A_3 C_1 X_4 + A_1 C_3 X_4)) \\
		
		D = (A_2 B_4 C_3 D1 - A_2 B_3 C_4 D1 - A_1 B_4 C_3 D_2 + A_1 B_3 C_4 D_2 - A_2 B_4 C_1 D_3 + A_1 B_4 C_2 D_3 + \\
A_2 B_1 C_4 D_3 - A_1 B_2 C_4 D_3 + A_4 (B_3 C_2 D1 - B_2 C_3 D1 - B_3 C_1 D_2 + B_1 C_3 D_2 + B_2 C_1 D_3 - B_1 C_2 D_3) + \\
A_2 B_3 C_1 D_4 - A_1 B_3 C_2 D_4 - A_2 B_1 C_3 D_4 + A_1 B_2 C_3 D_4 + A_3 (-B_4 C_2 D1 + B_2 C_4 D1 + B_4 C_1 D_2 - B_1 C_4 D_2 - B_2 C_1 D_4 + B_1 C_2 D_4))
		
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
