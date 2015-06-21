#include <stdlib.h>
#include <fstream>

#include "dynein_struct.h"


/* *********************************** DYNEIN FUNCTIONS ****************************************** */

void Dynein::set_state(states s) {
	state = s;
}

void Dynein::next_timestep() {
	
	r_blx = 0;
	r_mlx = 0;
	r_mrx = 0;
	r_brx = 0;

	r_bly = 0;
	r_mly = 0;
	r_mry = 0;
	r_bry = 0;

	f_blx = 0;
	f_mlx = 0;
	f_mrx = 0;
	f_brx = 0;

	f_bly = 0;
	f_mly = 0;
	f_mry = 0;
	f_bry = 0;
	
	double temp_d_bla;
	double temp_d_mla;
	double temp_d_mra;
	double temp_d_bra;
	
	if (state == LEFTBOUND) {
		temp_d_bla =
			(get_mly()*( (1/g)*(get_f_brx() + get_f_mlx() + get_f_mrx() + get_f_tx()) + get_r_brx() + get_r_mlx() + get_r_mrx() + get_r_tx() - get_d_mlx() - get_d_mrx() - get_d_tx()) -
			get_f_brx()*(1/g)*get_bly() - get_f_mlx()*(1/g)*get_bly() - get_f_mrx()*(1/g)*get_bly() - 
			get_f_tx()*(1/g)*get_bly() + get_f_bry()*(1/g)*(get_blx() - get_mlx()) + get_f_mly()*(1/g)*(get_blx() - get_mlx()) +
			get_f_mry()*(1/g)*get_blx() - get_f_mry()*(1/g)*get_mlx() + get_f_ty()*(1/g)*get_blx() - get_f_ty()*(1/g)*get_mlx() +
			ls*get_d_bra()*cos(bra)*(get_blx() - get_mlx()) + ls*get_d_bra()*sin(bra)*(get_bly() - get_mly()) +
			lt*(-(get_bly() - get_mly())*(get_d_mla()*sin(mla) - get_d_mra()*sin(mra)) + get_d_mla()*cos(mla)*(get_mlx() - get_blx()) + get_d_mra()*cos(mra)*(get_blx() - get_mlx())) -
			get_r_brx()*get_bly() - get_r_mlx()*get_bly() - get_r_mrx()*get_bly() - get_r_tx()*get_bly() +
			get_r_bry()*get_blx() - get_r_bry()*get_mlx() + get_r_mly()*get_blx() - get_r_mly()*get_mlx() +
			get_r_mry()*get_blx() - get_r_mry()*get_mlx() + get_r_ty()*get_blx() - get_r_ty()*get_mlx() - get_blx()*get_d_mly() -
			get_blx()*get_d_mry() - get_blx()*get_d_ty() + get_d_mlx()*get_bly() + get_d_mrx()*get_bly() + 
			get_d_tx()*get_bly() + get_mlx()*get_d_mly() + get_mlx()*get_d_mry() + get_mlx()*get_d_ty())
			/
			(ls*cos(bla)*(get_blx() - get_mlx()) + ls*sin(bla)*(get_bly() - get_mly()));
		
		temp_d_mla = 
			(get_ty()*((1/g)*(get_f_brx() + get_f_mrx() + get_f_tx()) + get_r_brx() + get_r_mrx() + get_r_tx() - get_d_mlx() - get_d_mrx() - get_d_tx()) -
			get_f_brx()*(1/g)*get_mly() - get_f_mrx()*(1/g)*get_mly() - get_f_tx()*(1/g)*get_mly() + get_f_bry()*(1/g)*(get_mlx() - get_tx()) + 
			get_f_mry()*(1/g)*(get_mlx() - get_tx()) + get_f_ty()*(1/g)*get_mlx() - get_f_ty()*(1/g)*get_tx() + 
			(get_mly() - get_ty()) * (ls*get_d_bra()*sin(bra) + lt*get_d_mra()*sin(mra)) +
			ls*get_d_bra()*cos(bra) * (get_mlx() - get_tx()) +
			lt*get_d_mra()*cos(mra) * (get_mlx() - get_tx()) -
			get_r_brx()*get_mly() - get_r_mrx()*get_mly() - get_r_tx()*get_mly() + get_r_bry()*get_mlx() -
			get_r_bry()*get_tx() + get_r_mry()*get_mlx() - get_r_mry()*get_tx() + get_r_ty()*get_mlx() -
			get_r_ty()*get_tx() + get_d_mlx()*get_mly() + get_d_mrx()*get_mly() + get_d_tx()*get_mly() - get_mlx()*get_d_mly() -
			get_mlx()*get_d_mry() -get_mlx()*get_d_ty() + get_tx()*get_d_mly() + get_tx()*get_d_mry() + get_tx()*get_d_ty())
			/ 
			(lt*cos(mla)*(get_mlx() - get_tx()) + lt*sin(mla)*(get_mly() - get_ty()));
			
		temp_d_mra =
			((get_mry() - get_ty())*((1/g)*(get_f_brx() + get_f_mrx()) + get_r_brx() + get_r_mrx() - get_d_mrx() - get_d_tx()) + 
			get_f_bry()*(1/g)*(get_tx() - get_mrx()) + get_f_mry()*(1/g)*(get_tx() - get_mrx()) + 
			ls*get_d_bra()*(cos(bra)*(get_tx() - get_mrx()) + sin(bra)*(get_ty() - get_mry())) -
			(get_mrx() - get_tx())*(get_r_bry() + get_r_bry() - get_d_mry() - get_d_ty()))
			/
			(lt*cos(mra)*(get_mrx()-get_tx()) + lt*sin(mra)*(get_mry()-get_ty()));
			
		temp_d_bra =
			((get_f_brx()*(1/g) + get_r_brx() - get_d_mrx())*(get_bry()-get_mry()) + (get_f_bry()*(1/g) + get_r_bry() - get_d_mry())*(get_mrx()-get_brx()))
			/
			(ls*cos(bra)*(get_brx()-get_mrx()) + ls*sin(bra)*(get_bry()-get_mry()));
		
		d_bla = temp_d_bla;
		d_mla = temp_d_mla;
		d_mra = temp_d_mra;
		d_bra = temp_d_bra;
		}
	else if (state == BOTHBOUND) {
		d_bla = 0;
		d_mla = 0;
		d_mra = 0;
		d_bra = 0;
	}
	else {
		d_bla = 0;
		d_mla = 0;
		d_mra = 0;
		d_bra = 0;
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


/*** Get coordinate velocities ***/

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
	if (state == LEFTBOUND) return lt * d_mra * sin(-mra) + get_d_tx();
	else if (state == BOTHBOUND) return -d_bra * ls * sin(bra);
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
	else if (state == BOTHBOUND) return d_bla * ls * cos(bla);
	else return lt * -d_mla * cos(-mla) + get_d_ty();
}

double Dynein::get_d_ty() {
	if (state == LEFTBOUND) return lt * d_mla * cos(mla) + get_d_mly();
	else if (state == BOTHBOUND) return lt/2 * (d_mra * cos(mra) + d_mla * cos(mla) + get_d_mly() + get_d_mry());
	else return lt * d_mra * cos(mra) + get_d_mry();
}

double Dynein::get_d_mry() {
	if (state == LEFTBOUND) return lt * -d_mra * cos(-mra) + get_d_ty();
	else if (state == BOTHBOUND) return d_bra * ls * cos(mra);
	else return ls * d_bra * cos(bra);
}

double Dynein::get_d_bry() {
	if (state == LEFTBOUND) return ls * -d_bra * cos(-bra) + get_d_mry();
	else if (state == BOTHBOUND) return 0;
	else return 0;
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
	return 0.1;
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


