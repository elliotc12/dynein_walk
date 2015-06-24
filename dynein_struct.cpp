#include <stdlib.h>
#include <fstream>

#include "dynein_struct.h"


/* *********************************** DYNEIN FUNCTIONS ****************************************** */

Dynein::Dynein(double bla_init, double mla_init, double mra_init, double bra_init) {
	blx = 0;
	bly = 0;
	
	d_bla = 0;
	d_mla = 0;
	d_mra = 0;
	d_bra = 0;
	
	d_blx = 0;    d_bly = 0;
	d_mlx = 0;    d_mly = 0;
	d_mrx = 0;    d_mry = 0;
	d_brx = 0;    d_bry = 0;
	
	bla = bla_init;
	mla = mla_init;
	mra = mra_init;
	bra = bra_init;
	
	//blx = blx;
	mlx = ls * cos(bla) + blx;
	tx  = lt * cos(mla) + mlx;
	mrx = lt * cos(-mra) + tx;
	brx = ls * cos(-bra) + mrx;
	
	//bly = bly;
	mly = ls * sin(bla) + bly;
	ty  = lt * sin(mla) + mly;
	mry = lt * sin(-mra) + ty;
	bry = ls * sin(-bra) + mry;
	
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
	
	double temp_d_bla;
	double temp_d_mla;
	double temp_d_mra;
	double temp_d_bra;
	
	double temp_blx;    double temp_bly;
	double temp_mlx;    double temp_mly;
	double temp_tx;     double temp_ty;
	double temp_mrx;    double temp_mry;
	double temp_brx;    double temp_bry;
	
	double temp_d_blx;    double temp_d_bly;
	double temp_d_mlx;    double temp_d_mly;
	double temp_d_tx;     double temp_d_ty;
	double temp_d_mrx;    double temp_d_mry;
	double temp_d_brx;    double temp_d_bry;
	
	if (state == LEFTBOUND) {
		
		temp_blx = blx;
		temp_mlx = ls * cos(bla) + temp_blx;
		temp_tx  = lt * cos(mla) + temp_mlx;
		temp_mrx = lt * cos(-mra) + temp_tx;
		temp_brx = ls * cos(-bra) + temp_mrx;
		
		temp_bly = bly;
	    temp_mly = ls * sin(bla) + temp_bly;
		temp_ty  = lt * sin(mla) + temp_mly;
		temp_mry = lt * sin(-mra) + temp_ty;
		temp_bry = ls * sin(-bra) + temp_mry;
		
		temp_d_bla =
			(mly*( (1/g)*(f_brx + f_mlx + f_mrx + f_tx) + r_brx + r_mlx + r_mrx + r_tx - d_mlx - d_mrx - d_tx) -
			f_brx*(1/g)*bly - f_mlx*(1/g)*bly - f_mrx*(1/g)*bly - f_tx*(1/g)*bly + f_bry*(1/g)*(blx - mlx) + f_mly*(1/g)*(blx - mlx) +
			f_mry*(1/g)*blx - f_mry*(1/g)*mlx + f_ty*(1/g)*blx - f_ty*(1/g)*mlx + ls*d_bra*cos(bra)*(blx - mlx) + ls*d_bra*sin(bra)*(bly - mly) +
			lt*(-(bly - mly)*(d_mla*sin(mla) - d_mra*sin(mra)) + d_mla*cos(mla)*(mlx - blx) + d_mra*cos(mra)*(blx - mlx)) -
			r_brx*bly - r_mlx*bly - r_mrx*bly - r_tx*bly + r_bry*blx - r_bry*mlx + r_mly*blx - r_mly*mlx + r_mry*blx - r_mry*mlx + r_ty*blx -
			r_ty*mlx - blx*d_mly - blx*d_mry - blx*d_ty + d_mlx*bly + d_mrx*bly + d_tx*bly + mlx*d_mly + mlx*d_mry + mlx*d_ty)
			/
			(ls*cos(bla)*(blx - mlx) + ls*sin(bla)*(bly - mly));

		temp_d_mla = 
			(ty*((1/g)*(f_brx + f_mrx + f_tx) + r_brx + r_mrx + r_tx - d_mlx - d_mrx - d_tx) - f_brx*(1/g)*mly - f_mrx*(1/g)*mly - f_tx*(1/g)*mly + 
			f_bry*(1/g)*(mlx - tx) + f_mry*(1/g)*(mlx - tx) + f_ty*(1/g)*mlx - f_ty*(1/g)*tx + (mly - ty) * (ls*d_bra*sin(bra) + 
			lt*d_mra*sin(mra)) + ls*d_bra*cos(bra) * (mlx - tx) + lt*d_mra*cos(mra) * (mlx - tx) - r_brx*mly - r_mrx*mly - r_tx*mly + r_bry*mlx - r_bry*tx + 
			r_mry*mlx - r_mry*tx + r_ty*mlx - r_ty*tx + d_mlx*mly + d_mrx*mly + d_tx*mly - mlx*d_mly - mlx*d_mry -mlx*d_ty + tx*d_mly + tx*d_mry + tx*d_ty)
			/ 
			(lt*cos(mla)*(mlx - tx) + lt*sin(mla)*(mly - ty));

		temp_d_mra =
			((mry - ty)*((1/g)*(f_brx + f_mrx) + r_brx + r_mrx - d_mrx - d_tx) + f_bry*(1/g)*(tx - mrx) + f_mry*(1/g)*(tx - mrx) + 
			ls*d_bra*(cos(bra)*(tx - mrx) + sin(bra)*(ty - mry)) - (mrx - tx)*(r_bry + r_bry - d_mry - d_ty))
			/
			(lt*cos(mra)*(mrx-tx) + lt*sin(mra)*(mry-ty));

		temp_d_bra =
			((f_brx*(1/g) + r_brx - d_mrx)*(bry-mry) + (f_bry*(1/g) + r_bry - d_mry)*(mrx-brx))
			/
			(ls*cos(bra)*(brx-mrx) + ls*sin(bra)*(bry-mry));
		
		d_bla = temp_d_bla;
		d_mla = temp_d_mla;
		d_mra = temp_d_mra;
		d_bra = temp_d_bra;
		
		temp_d_blx = 0;
		temp_d_mlx = ls * d_bla * -sin(bla);
		temp_d_tx  = lt * d_mla * -sin(mla) + temp_d_mlx;
		temp_d_mrx = lt * d_mra * sin(-mra) + temp_d_tx;
		temp_d_brx = ls * d_bra * sin(-bra) + temp_d_mrx;
		
		temp_d_bly = 0;
		temp_d_mly = ls * d_bla * cos(bla);
		temp_d_ty  = lt * d_mla * cos(mla) + temp_d_mly;
		temp_d_mry = lt * -d_mra * cos(-mra) + temp_d_ty;
		temp_d_bry = ls * -d_bra * cos(-bra) + temp_d_mry;
		
		d_blx = temp_d_blx;      d_bly = temp_d_bly;
		d_mlx = temp_d_mlx;      d_mly = temp_d_mly;
		d_tx  = temp_d_tx;       d_ty  = temp_d_ty; 
		d_mrx = temp_d_mrx;      d_mry = temp_d_mry;
		d_brx = temp_d_brx;      d_bry = temp_d_bry;
		
		blx = temp_blx;      bly = temp_bly;
		mlx = temp_mlx;      mly = temp_mly;
		tx  = temp_tx;       ty  = temp_ty; 
		mrx = temp_mrx;      mry = temp_mry;
		brx = temp_brx;      bry = temp_bry;
		
		//printf("d_bla: %f\n", );
		
	} else if (state == BOTHBOUND) {
		
		temp_d_blx = 0;
		temp_d_mlx = -d_bla * ls * sin(bla);
		temp_d_mrx = -d_bra * ls * sin(bra);
		temp_d_tx  = lt/2 * (-d_mra * sin(mra) + -d_mla * sin(mla) + temp_d_mlx + temp_d_mrx);
		temp_d_brx = 0;
		            
	    temp_d_bly = 0;
		temp_d_mly = d_bla * ls * cos(bla);
		temp_d_mry = d_bra * ls * cos(mra);
		temp_d_ty  = lt/2 * (d_mra * cos(mra) + d_mla * cos(mla) + temp_d_mly + temp_d_mry);
		temp_d_bry = 0;
		
		d_bla = 0;
		d_mla = 0;
		d_mra = 0;
		d_bra = 0;
		
		temp_blx = 0;      temp_bly = 0;
		temp_mlx = 0;      temp_mly = 0;
		temp_tx  = 0;      temp_ty  = 0;
		temp_mrx = 0;      temp_mry = 0;
		temp_brx = 0;      temp_bry = 0;
		
		blx = temp_blx;      bly = temp_bly;
		mlx = temp_mlx;      mly = temp_mly;
		tx  = temp_tx;       ty  = temp_ty; 
		mrx = temp_mrx;      mry = temp_mry;
		brx = temp_brx;      bry = temp_bry;
	} else {
		
		temp_d_brx = 0;
		temp_d_mrx = ls * d_bra * -sin(bra);
		temp_d_tx  = lt * d_mra * -sin(mra) + temp_d_mrx;
		temp_d_mlx = lt * d_mla * sin(-mla) + temp_d_tx;
		temp_d_blx = ls * d_bla * sin(-bla) + temp_d_mlx;
		
		temp_d_bry = 0;
		temp_d_mry = ls * d_bra * cos(bra);
		temp_d_ty  = lt * d_mra * cos(mra) + temp_d_mry;
		temp_d_mly = lt * -d_mla * cos(-mla) + temp_d_ty;
		temp_d_bly = ls * -d_bla * cos(-bla) + temp_d_mly;
		
		d_bla = 0;
		d_mla = 0;
		d_mra = 0;
		d_bra = 0;
		
		temp_blx = 0;      temp_bly = 0;
		temp_mlx = 0;      temp_mly = 0;
		temp_tx  = 0;      temp_ty  = 0;
		temp_mrx = 0;      temp_mry = 0;
		temp_brx = 0;      temp_bry = 0;
		
		blx = temp_blx;      bly = temp_bly;
		mlx = temp_mlx;      mly = temp_mly;
		tx  = temp_tx;       ty  = temp_ty; 
		mrx = temp_mrx;      mry = temp_mry;
		brx = temp_brx;      bry = temp_bry;
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
	fprintf(data_file, "%.6f\t%12.6f\t%12.6f\t%.3f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%d\n", 
	get_KE(), get_PE(), get_KE() + get_PE(), t, blx, bly, mlx, mly, tx, ty, mrx, mry, brx, bry, get_state());
	fclose(data_file);
}
