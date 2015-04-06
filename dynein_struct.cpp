#include <stdlib.h>
#include <fstream>

#include "dynein_struct.h"


/* *********************************** DYNEIN FUNCTIONS ****************************************** */

void Dynein::set_state(states s) {
	state = s;
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



/*** Get angular velocities ***/

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



/*** Get coordinate velocities ***/

double Dynein::get_d_blx() {
	if (state == LEFTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
	else return ls * d_bla * sin(-bla) + get_d_mlx();
}

double Dynein::get_d_mlx() {
	if (state == LEFTBOUND) return ls * d_bla * -sin(bla);
	else if (state == BOTHBOUND) return 0;
	else return lt * d_mla * sin(-mla) + get_d_tx();
}

double Dynein::get_d_tx() {
	if (state == LEFTBOUND) return lt * d_mla * -sin(mla) + get_d_mlx();
	else if (state == BOTHBOUND) return 0;
	else return lt * d_mra * -sin(mra) + get_d_mrx();
}

double Dynein::get_d_mrx() {
	if (state == LEFTBOUND) return lt * d_mra * sin(-mra) + get_d_tx();
	else if (state == BOTHBOUND) return 0;
	else return ls * d_bra * -sin(bra);
}

double Dynein::get_d_brx() {
	if (state == LEFTBOUND) return ls * d_bra * sin(-bra) + get_d_mrx();
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
	else if (state == BOTHBOUND) return 0;
	else return lt * -d_mla * cos(-mla) + get_d_ty();
}

double Dynein::get_d_ty() {
	if (state == LEFTBOUND) return lt * d_mla * cos(mla) + get_d_mly();
	else if (state == BOTHBOUND) return 0;
	else return lt * d_mra * cos(mra) + get_d_mry();
}

double Dynein::get_d_mry() {
	if (state == LEFTBOUND) return lt * -d_mra * cos(-mra) + get_d_ty();
	else if (state == BOTHBOUND) return 0;
	else return ls * d_bra * cos(bra);
}

double Dynein::get_d_bry() {
	if (state == LEFTBOUND) return ls * -d_bra * cos(-bra) + get_d_mry();
	else if (state == BOTHBOUND) return 0;
	else return 0;
}



/*** Get cartesian accelerations ***/

double Dynein::get_dd_blx() {
	return 0;
}

double Dynein::get_dd_mlx() {
	return 0;
}

double Dynein::get_dd_mrx() {
	return 0;
}

double Dynein::get_dd_brx() {
	return 0;
}

double Dynein::get_dd_bly() {
	return 0;
}

double Dynein::get_dd_mly() {
	return 0;
}

double Dynein::get_dd_mry() {
	return 0;
}

double Dynein::get_dd_bry() {
	return 0;
}



/*** Get coordinates ***/

double Dynein::get_blx() {
	if (state == LEFTBOUND) return blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_bly(){
	if (state == LEFTBOUND) return bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mlx() {
	if (state == LEFTBOUND) return ls * cos(get_bla()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mly(){
	if (state == LEFTBOUND) return ls * sin(get_bla()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_tx() {
	if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_ty(){
	if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mrx() {
	if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + lt * cos(-get_mra()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mry(){
	if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + lt * sin(-get_mra()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_brx() {
	if (state == LEFTBOUND) return ls * cos(get_bla()) + lt * cos(get_mla()) + lt * cos(-get_mra()) + ls * cos(-get_bra()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_bry(){
	if (state == LEFTBOUND) return ls * sin(get_bla()) + lt * sin(get_mla()) + lt * sin(-get_mra()) + ls * sin(-get_bra()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
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
	return (1.0/2) * kbl * pow(get_bla() - ba, 2) + (1.0/2) * kml * pow(M_PI - get_bla() + get_mla() - ma, 2) + (1.0/2) * kt * 
		pow(M_PI - get_mla() - get_mra() - ta, 2) + (1.0/2) * kmr * pow(M_PI - get_bra() + get_mra() - ma, 2);
}

double Dynein::get_KE() {
	return 1.0/2.0 * mb * (square(- (ls * get_d_bla() * sin(get_bla())) - ls * get_d_bra() * sin(get_bra()) - lt * get_d_mla() * sin(get_mla()) 
		- lt * get_d_mra() * sin(get_mra())) +  square(ls * get_d_bla() * cos(get_bla()) -  ls * get_d_bra() * cos(get_bra()) + lt * get_d_mla() 
		* cos(get_mla()) - lt * get_d_mra() * cos(get_mra()))) + 1.0/2.0 * mm * (square(ls) * square(get_d_bla()) * square(sin(get_bla())) 
		+ square(ls) * square(get_d_bla()) * square(cos(get_bla()))) + 1.0/2.0 * mm * (square(- (ls * get_d_bla() * sin(get_bla())) - lt * 
		get_d_mla() * sin(get_mla()) - lt * get_d_mra() * sin(get_mra())) +  square(ls * get_d_bla() * cos(get_bla()) +  lt * get_d_mla() * 
		cos(get_mla()) - lt * get_d_mra() * cos(get_mra()))) + 1.0/2.0 * mt * (square(- (ls * get_d_bla() * sin(get_bla())) - lt * get_d_mla() 
		* sin(get_mla())) +  square(ls * get_d_bla() * cos(get_bla()) +  lt * get_d_mla() * cos(get_mla())));
	
	//Rewrite with new d_blx, etc.
}

void Dynein::log(double t) {
	FILE* data_file = fopen("data.txt", "a+");
	fprintf(data_file, "%.6f\t%12.6f\t%12.6f\t%.3f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%d\n", 
	get_KE(), get_PE(), get_KE() + get_PE(), t, get_blx(), get_bly(), get_mlx(), get_mly(), get_tx(), get_ty(), get_mrx(), 
	get_mry(), get_brx(), get_bry(), get_state());
	fclose(data_file);
}


