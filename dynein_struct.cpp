#include <stdlib.h>
#include <fstream>

#include "dynein_struct.h"

double fblx = 0;	double fbly = 0;
double fmlx = 0;	double fmly = 0;
double ftx = 0;		double fty = 0;
double fmrx = 0;	double fmry = 0;
double fbrx = 0;	double fbry = 0;


/* *********************************** DYNEIN FUNCTIONS ****************************************** */

/*** Accessor Methods ***/

void Dynein::set_state(states s) {
	state = s;
}

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

double Dynein::get_d_blx() {
	return d_bla;
}

double Dynein::get_d_mlx() {
	return d_mla;
}

double Dynein::get_d_mrx() {
	return d_mra;
}

double Dynein::get_d_brx() {
	return d_bra;
}

double Dynein::get_d_bly() {
	return d_bla;
}

double Dynein::get_d_mly() {
	return d_mla;
}

double Dynein::get_d_mry() {
	return d_mra;
}

double Dynein::get_d_bry() {
	return d_bra;
}

double Dynein::get_dd_blx() {
	return d_bla;
}

double Dynein::get_dd_mlx() {
	return d_mla;
}

double Dynein::get_dd_mrx() {
	return d_mra;
}

double Dynein::get_dd_brx() {
	return d_bra;
}

double Dynein::get_dd_bly() {
	return d_bla;
}

double Dynein::get_dd_mly() {
	return d_mla;
}

double Dynein::get_dd_mry() {
	return d_mra;
}

double Dynein::get_dd_bry() {
	return d_bra;
}

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





/* *********************************** UTILITY FUNCTIONS ****************************************** */

double dist(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1) * range;
}

void resetLog(Dynein* dyn) {
	FILE* data_file = fopen("data.txt", "w");
	FILE* config_file = fopen("config.txt", "w");
	
	fprintf(config_file, "#inctime\truntime\tstate\n%+.3f\t%+.3f\t%d\n", inctime, runtime, (int) dyn->get_state());
	fprintf(data_file,
		"#KE\t\t\t\tPE\t\t\t\tEnergy\t\tt\t\tblX\t\t\tblY\t\t\tmlX\t\t\tmlY\t\t\ttX\t\t\ttY\t\t\tmrX\t\t\tmrY\t\t\tbrX\t\t\tbrY\t\t\tS\n");
	
	fclose(data_file);
	fclose(config_file);
}

double square(double num) {
	return num * num;
}

double cube(double num) {
	return num * num * num;
}

double fifth(double num) {
	return num * num * num * num * num;
}
