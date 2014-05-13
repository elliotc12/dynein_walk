//Do testing in Mathematica to numerically evaluate each acceleration equation and compare those results to those of the replace.py'd C a. functions
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

const double inctime = 0.1;
const double runtime = 1000.0;

const double Lt = 10.0;
const double Ls = 10.0;

const double kt  = 1.0; // Higher spring constant -> less deviation from equilibrium
const double kml = 1.0;
const double kmr = 1.0;
const double kbl = 1.0;
const double kbr = 1.0;
	
const double mb = 1.0; // Higher mass constant -> less movement
const double mm = 1.0;
const double mt = 1.0;

const double bA = (108.0 / 180) * M_PI;
const double mA = (108.0 / 180) * M_PI;
const double tA = (108.0 / 180) * M_PI;

const double bla_init = (108.0 / 180) * M_PI;
const double mla_init = (36.0 / 180) * M_PI;
const double mra_init = (36.0 / 180) * M_PI;
const double bra_init = (100.0 / 180) * M_PI;

FILE* data_file;
FILE* config_file;

typedef enum {
	LEFTBOUND,
	RIGHTBOUND,
	BOTHBOUND
} states;

class Dynein {
	public:
		void set_bla(double d);
		void set_mla(double d);
		void set_mra(double d);
		void set_bra(double d);
		
		void set_blx(double d);
		void set_bly(double d);
		
		void set_state(states s);
		
		void set_d_bla(double d);
		void set_d_mla(double d);
		void set_d_mra(double d);
		void set_d_bra(double d);
		
		double get_bly();
		double get_blx();
		double get_mlx();
		double get_mly();
		double get_tx();
		double get_ty();
		double get_mrx();
		double get_mry();
		double get_brx();
		double get_bry();
		
		double get_bla();
		double get_mla();
		double get_mra();
		double get_bra();
		
		double get_d_bla();
		double get_d_mla();
		double get_d_mra();
		double get_d_bra();
		
		double get_dd_bla();
		double get_dd_mla();
		double get_dd_mra();
		double get_dd_bra();
		
		double get_PE();
		double get_KE();
		
		states get_state();
		
		void log(double t);
		
	private:
		double bla;
		double mla;
		double mra;
		double bra;
	
		double d_bla; //Angular Velocities
		double d_mla;
		double d_mra;
		double d_bra;
		
		double blx;
		double bly;
		
		states state;
};

double randAngle(double range);
double dist(double d, double h, double i, double j);
void   resetLog(Dynein* dyn);
double square(double num);
double cube(double num);

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
	if (state == LEFTBOUND) return Ls * cos(get_bla()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mly(){
	if (state == LEFTBOUND) return Ls * sin(get_bla()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_tx() {
	if (state == LEFTBOUND) return Ls * cos(get_bla()) + Lt * cos(get_mla()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_ty(){
	if (state == LEFTBOUND) return Ls * sin(get_bla()) + Lt * sin(get_mla()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mrx() {
	if (state == LEFTBOUND) return Ls * cos(get_bla()) + Lt * cos(get_mla()) + Lt * cos(-get_mra()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_mry(){
	if (state == LEFTBOUND) return Ls * sin(get_bla()) + Lt * sin(get_mla()) + Lt * sin(-get_mra()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_brx() {
	if (state == LEFTBOUND) return Ls * cos(get_bla()) + Lt * cos(get_mla()) + Lt * cos(-get_mra()) + Ls * cos(-get_bra()) + blx;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_bry(){
	if (state == LEFTBOUND) return Ls * sin(get_bla()) + Lt * sin(get_mla()) + Lt * sin(-get_mra()) + Ls * sin(-get_bra()) + bly;
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_dd_bla() { 
	if (state == LEFTBOUND) return 
		- (4.0 * square(Ls) * Lt * mt * sin(2.0 * (get_bla() -  get_mla())) * square(get_d_bla()) * square(mm) + 8.0 * Ls * square(Lt) * mt 
		* sin(get_bla() -  get_mla()) * square(get_d_mla()) * square(mm) + 4.0 * Ls * square(Lt) * mt * sin(get_bla() -  2.0 * get_mla() - 
		 get_mra()) * square(get_d_mra()) * square(mm) + 4.0 * Ls * square(Lt) * mt * sin(get_bla() +  get_mra()) * square(get_d_mra()) * square(mm) 
		-  4.0 * bA * kbl * Lt * square(mm) +  4.0 * kml * Lt * mA * square(mm) + 4.0 * kml * Ls * mA * cos(get_bla() -  get_mla()) * square(mm) 
		- 4.0 * kt * Ls * tA * cos(get_bla() -  get_mla()) * square(mm) - 4.0 * kml * Ls * M_PI * cos(get_bla() -  get_mla()) * square(mm) 
		+ 4.0 * kt * Ls * M_PI * cos(get_bla() -  get_mla()) * square(mm) + 4.0 * kmr * Ls * mA * cos(get_bla() -  2.0 * get_mla() -  get_mra()) 
		* square(mm) - 4.0 * kt * Ls * tA * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * square(mm) - 4.0 * kmr * Ls * M_PI * cos(get_bla() 
		-  2.0 * get_mla() -  get_mra()) * square(mm) + 4.0 * kt * Ls * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * square(mm) 
		- 4.0 * kmr * Ls * mA * cos(get_bla() +  get_mra()) * square(mm) + 4.0 * kt * Ls * tA * cos(get_bla() +  get_mra()) * square(mm) + 
		4.0 * kmr * Ls * M_PI * cos(get_bla() +  get_mra()) * square(mm) - 4.0 * kt * Ls * M_PI * cos(get_bla() +  get_mra()) * square(mm) 
		+ 4.0 * bA * kbl * Lt * cos(2.0 * (get_mla() +  get_mra())) * square(mm) - 4.0 * kml * Lt * mA * cos(2.0 * (get_mla() +  get_mra())) 
		* square(mm) + 4.0 * kml * Lt * M_PI * cos(2.0 * (get_mla() +  get_mra())) * square(mm) - 4.0 * kml * Ls * mA * cos(get_bla() +  get_mla() 
		+  2.0 * get_mra()) * square(mm) + 4.0 * kt * Ls * tA * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) + 4.0 * kml * Ls 
		* M_PI * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) - 4.0 * kt * Ls * M_PI * cos(get_bla() +  get_mla() +  2.0 * get_mra()) 
		* square(mm) - 4.0 * kt * Ls * cos(get_bla() -  get_mla()) * get_mra() * square(mm) - 4.0 * kmr * Ls * cos(get_bla() -  2.0 * get_mla() 
		-  get_mra()) * get_mra() * square(mm) - 4.0 * kt * Ls * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mra() * square(mm) + 
		4.0 * kmr * Ls * cos(get_bla() +  get_mra()) * get_mra() * square(mm) + 4.0 * kt * Ls * cos(get_bla() +  get_mra()) * get_mra() * square(mm) 
		+ 4.0 * kt * Ls * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * get_mra() * square(mm) + 8.0 * kmr * Ls * get_bra() * sin(get_bla() 
		-  get_mla()) * sin(get_mla() +  get_mra()) * square(mm) -  4.0 * kml * Lt * M_PI * square(mm) + 4.0 * square(Ls) * Lt * square(mt) 
		* sin(2.0 * (get_bla() -  get_mla())) * square(get_d_bla()) * mm + 4.0 * square(Ls) * Lt * mb * mt * sin(2.0 * (get_bla() -  get_mla())) 
		* square(get_d_bla()) * mm + 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() +  get_bra()) * square(get_d_bra()) * mm + 2.0 * square(Ls) 
		* Lt * mb * mt * sin(get_bla() -  get_bra() -  2.0 * get_mla()) * square(get_d_bra()) * mm + 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() 
		-  get_bra() +  2.0 * get_mra()) * square(get_d_bra()) * mm + 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() +  get_bra() -  2.0 * 
		(get_mla() +  get_mra())) * square(get_d_bra()) * mm + 8.0 * Ls * square(Lt) * square(mt) * sin(get_bla() -  get_mla()) * square(get_d_mla()) 
		* mm + 8.0 * Ls * square(Lt) * mb * mt * sin(get_bla() -  get_mla()) * square(get_d_mla()) * mm + 4.0 * Ls * square(Lt) * mb * mt * 
		sin(get_bla() -  2.0 * get_mla() -  get_mra()) * square(get_d_mra()) * mm + 4.0 * Ls * square(Lt) * mb * mt * sin(get_bla() +  get_mra()) 
		* square(get_d_mra()) * mm -  4.0 * bA * kbl * Lt * mb * mm +  4.0 * kml * Lt * mA * mb * mm - 8.0 * bA * kbl * Lt * mt * mm +  8.0 
		* kml * Lt * mA * mt * mm + 2.0 * kmr * Lt * mA * mt * cos(get_bla() +  get_bra()) * mm - 2.0 * kmr * Lt * mt * M_PI * cos(get_bla() 
		+  get_bra()) * mm - 2.0 * kmr * Lt * mA * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) * mm + 2.0 * kmr * Lt * mt * M_PI * 
		cos(get_bla() -  get_bra() -  2.0 * get_mla()) * mm + 4.0 * kml * Ls * mA * mb * cos(get_bla() -  get_mla()) * mm + 8.0 * kml * Ls 
		* mA * mt * cos(get_bla() -  get_mla()) * mm - 4.0 * kt * Ls * mb * tA * cos(get_bla() -  get_mla()) * mm - 8.0 * kt * Ls * mt * tA 
		* cos(get_bla() -  get_mla()) * mm - 4.0 * kml * Ls * mb * M_PI * cos(get_bla() -  get_mla()) * mm + 4.0 * kt * Ls * mb * M_PI * cos(get_bla() 
		-  get_mla()) * mm - 8.0 * kml * Ls * mt * M_PI * cos(get_bla() -  get_mla()) * mm + 8.0 * kt * Ls * mt * M_PI * cos(get_bla() -  get_mla()) 
		* mm + 4.0 * kmr * Ls * mA * mb * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm + 4.0 * kmr * Ls * mA * mt * cos(get_bla() - 
		 2.0 * get_mla() -  get_mra()) * mm - 4.0 * kt * Ls * mb * tA * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm - 4.0 * kt * Ls 
		* mt * tA * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm - 4.0 * kmr * Ls * mb * M_PI * cos(get_bla() -  2.0 * get_mla() -  
		get_mra()) * mm + 4.0 * kt * Ls * mb * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm - 4.0 * kmr * Ls * mt * M_PI * cos(get_bla() 
		-  2.0 * get_mla() -  get_mra()) * mm + 4.0 * kt * Ls * mt * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm - 4.0 * kmr 
		* Ls * mA * mb * cos(get_bla() +  get_mra()) * mm - 4.0 * kmr * Ls * mA * mt * cos(get_bla() +  get_mra()) * mm + 4.0 * kt * Ls * mb 
		* tA * cos(get_bla() +  get_mra()) * mm + 4.0 * kt * Ls * mt * tA * cos(get_bla() +  get_mra()) * mm + 4.0 * kmr * Ls * mb * M_PI * 
		cos(get_bla() +  get_mra()) * mm - 4.0 * kt * Ls * mb * M_PI * cos(get_bla() +  get_mra()) * mm + 4.0 * kmr * Ls * mt * M_PI * cos(get_bla() 
		+  get_mra()) * mm - 4.0 * kt * Ls * mt * M_PI * cos(get_bla() +  get_mra()) * mm + 4.0 * bA * kbl * Lt * mb * cos(2.0 * (get_mla() 
		+  get_mra())) * mm - 4.0 * kml * Lt * mA * mb * cos(2.0 * (get_mla() +  get_mra())) * mm + 4.0 * kml * Lt * mb * M_PI * cos(2.0 * 
		(get_mla() +  get_mra())) * mm - 2.0 * kmr * Lt * mA * mt * cos(get_bla() -  get_bra() +  2.0 * get_mra()) * mm + 2.0 * kmr * Lt * 
		mt * M_PI * cos(get_bla() -  get_bra() +  2.0 * get_mra()) * mm - 4.0 * kml * Ls * mA * mb * cos(get_bla() +  get_mla() +  2.0 * get_mra()) 
		* mm + 4.0 * kt * Ls * mb * tA * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * mm + 4.0 * kml * Ls * mb * M_PI * cos(get_bla() + 
		 get_mla() +  2.0 * get_mra()) * mm - 4.0 * kt * Ls * mb * M_PI * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * mm + 2.0 * kmr * 
		Lt * mA * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * mm - 2.0 * kmr * Lt * mt * M_PI * cos(get_bla() +  get_bra() 
		-  2.0 * (get_mla() +  get_mra())) * mm - 2.0 * kmr * Lt * mt * cos(get_bla() +  get_bra()) * get_mra() * mm + 2.0 * kmr * Lt * mt 
		* cos(get_bla() -  get_bra() -  2.0 * get_mla()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(get_bla() -  get_mla()) * get_mra() * 
		mm - 8.0 * kt * Ls * mt * cos(get_bla() -  get_mla()) * get_mra() * mm - 4.0 * kmr * Ls * mb * cos(get_bla() -  2.0 * get_mla() -  
		get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mra() * mm - 4.0 * kmr * Ls 
		* mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) 
		* get_mra() * mm + 4.0 * kmr * Ls * mb * cos(get_bla() +  get_mra()) * get_mra() * mm + 4.0 * kt * Ls * mb * cos(get_bla() +  get_mra()) 
		* get_mra() * mm + 4.0 * kmr * Ls * mt * cos(get_bla() +  get_mra()) * get_mra() * mm + 4.0 * kt * Ls * mt * cos(get_bla() +  get_mra()) 
		* get_mra() * mm + 2.0 * kmr * Lt * mt * cos(get_bla() -  get_bra() +  2.0 * get_mra()) * get_mra() * mm + 4.0 * kt * Ls * mb * cos(get_bla() 
		+  get_mla() +  2.0 * get_mra()) * get_mra() * mm - 2.0 * kmr * Lt * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) 
		* get_mra() * mm -  4.0 * kmr * Lt * mt * get_bra() * sin(get_bla() -  get_mla()) * sin(get_bra() +  get_mla()) * mm - 4.0 * kmr * 
		Lt * mt * get_bra() * sin(get_bla() -  get_mla()) * sin(get_bra() -  get_mla() -  2.0 * get_mra()) * mm + 8.0 * kmr * Ls * mb * get_bra() 
		* sin(get_bla() -  get_mla()) * sin(get_mla() +  get_mra()) * mm + 8.0 * kmr * Ls * mt * get_bra() * sin(get_bla() -  get_mla()) * 
		sin(get_mla() +  get_mra()) * mm -  4.0 * kml * Lt * mb * M_PI * mm - 8.0 * kml * Lt * mt * M_PI * mm + 2.0 * square(Ls) * Lt * mb 
		* square(mt) * sin(2.0 * (get_bla() -  get_mla())) * square(get_d_bla()) -  square(Ls) * Lt * mb * square(mt) * sin(2.0 * (get_bla() 
		+  get_bra() -  get_mla() -  get_mra())) * square(get_d_bla()) -  square(Ls) * Lt * mb * square(mt) * sin(2.0 * (get_bla() -  get_bra() 
		- get_mla() +  get_mra())) * square(get_d_bla()) +  4.0 * Ls * square(Lt) * mb * square(mt) * sin(get_bla() -  get_mla()) * square(get_d_mla()) 
		- 2.0 * Ls * square(Lt) * mb * square(mt) * sin(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) * square(get_d_mla()) 
		-  2.0 * Ls * square(Lt) * mb * square(mt) * sin(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) * square(get_d_mla()) 
		-  4.0 * bA * kbl * Lt * mb * mt +  4.0 * kml * Lt * mA * mb * mt + 2.0 * kmr * Lt * mA * mb * mt * cos(get_bla() +  get_bra()) - 2.0 
		* kmr * Lt * mb * mt * M_PI * cos(get_bla() +  get_bra()) - 2.0 * kmr * Lt * mA * mb * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) 
		+ 2.0 * kmr * Lt * mb * mt * M_PI * cos(get_bla() -  get_bra() -  2.0 * get_mla()) + 4.0 * kml * Ls * mA * mb * mt * cos(get_bla() 
		-  get_mla()) - 4.0 * kt * Ls * mb * mt * tA * cos(get_bla() -  get_mla()) - 4.0 * kml * Ls * mb * mt * M_PI * cos(get_bla() -  get_mla()) 
		+ 4.0 * kt * Ls * mb * mt * M_PI * cos(get_bla() -  get_mla()) - 2.0 * kml * Ls * mA * mb * mt * cos(get_bla() +  2.0 * get_bra() - 
		 get_mla() -  2.0 * get_mra()) + 2.0 * kt * Ls * mb * mt * tA * cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) + 
		2.0 * kml * Ls * mb * mt * M_PI * cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) - 2.0 * kt * Ls * mb * mt * M_PI 
		* cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) + 4.0 * bA * kbl * Lt * mb * mt * cos(2.0 * (get_bra() -  get_mra())) 
		- 4.0 * kml * Lt * mA * mb * mt * cos(2.0 * (get_bra() -  get_mra())) + 4.0 * kml * Lt * mb * mt * M_PI * cos(2.0 * (get_bra() -  get_mra())) 
		+ 2.0 * kmr * Ls * mA * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mra()) - 2.0 * kt * Ls * mb * mt * tA * cos(get_bla() +  
		2.0 * get_bra() -  get_mra()) - 2.0 * kmr * Ls * mb * mt * M_PI * cos(get_bla() +  2.0 * get_bra() -  get_mra()) + 2.0 * kt * Ls * 
		mb * mt * M_PI * cos(get_bla() +  2.0 * get_bra() -  get_mra()) + 2.0 * kmr * Ls * mA * mb * mt * cos(get_bla() -  2.0 * get_mla() 
		- get_mra()) - 2.0 * kt * Ls * mb * mt * tA * cos(get_bla() -  2.0 * get_mla() -  get_mra()) - 2.0 * kmr * Ls * mb * mt * M_PI * cos(get_bla() 
		-  2.0 * get_mla() -  get_mra()) + 2.0 * kt * Ls * mb * mt * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) - 2.0 * kmr * Ls 
		* mA * mb * mt * cos(get_bla() +  get_mra()) + 2.0 * kt * Ls * mb * mt * tA * cos(get_bla() +  get_mra()) + 2.0 * kmr * Ls * mb * mt 
		* M_PI * cos(get_bla() +  get_mra()) - 2.0 * kt * Ls * mb * mt * M_PI * cos(get_bla() +  get_mra()) - 2.0 * kmr * Ls * mA * mb * mt 
		* cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) + 2.0 * kt * Ls * mb * mt * tA * cos(get_bla() -  2.0 * get_bra() 
		-  2.0 * get_mla() +  get_mra()) + 2.0 * kmr * Ls * mb * mt * M_PI * cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) 
		- 2.0 * kt * Ls * mb * mt * M_PI * cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) - 2.0 * kmr * Lt * mA * mb * mt 
		* cos(get_bla() -  get_bra() +  2.0 * get_mra()) + 2.0 * kmr * Lt * mb * mt * M_PI * cos(get_bla() -  get_bra() +  2.0 * get_mra()) 
		- 2.0 * kml * Ls * mA * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) + 2.0 * kt * Ls * mb * mt * tA 
		* cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) + 2.0 * kml * Ls * mb * mt * M_PI * cos(get_bla() -  2.0 * get_bra() 
		-  get_mla() +  2.0 * get_mra()) - 2.0 * kt * Ls * mb * mt * M_PI * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) 
		+ get_bla() * (4.0 * kml * Ls * (mb * (mm +  mt) +  mm * (mm +  2.0 * mt)) * cos(get_bla() -  get_mla()) - 2.0 * (- 2.0 * kbl * Lt 
		* square(mm) -  2.0 * kml * Lt * square(mm) + 2.0 * kbl * Lt * cos(2.0 * (get_mla() +  get_mra())) * square(mm) + 2.0 * kml * Lt * 
		cos(2.0 * (get_mla() +  get_mra())) * square(mm) + 2.0 * kml * Ls * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) - 2.0 
		* kbl * Lt * mb * mm -  2.0 * kml * Lt * mb * mm -  4.0 * kbl * Lt * mt * mm - 4.0 * kml * Lt * mt * mm + 2.0 * kbl * Lt * mb * cos(2.0 
		* (get_mla() +  get_mra())) * mm + 2.0 * kml * Lt * mb * cos(2.0 * (get_mla() +  get_mra())) * mm + 2.0 * kml * Ls * mb * cos(get_bla() 
		+  get_mla() +  2.0 * get_mra()) * mm - 2.0 * kbl * Lt * mb * mt -  2.0 * kml * Lt * mb * mt + kml * Ls * mb * mt * cos(get_bla() + 
		 2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) + 2.0 * (kbl +  kml) * Lt * mb * mt * cos(2.0 * (get_bra() -  get_mra())) + kml * 
		Ls * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()))) + 2.0 * kmr * Lt * mA * mb * mt * cos(get_bla() 
		+  get_bra() -  2.0 * (get_mla() +  get_mra())) - 2.0 * kmr * Lt * mb * mt * M_PI * cos(get_bla() +  get_bra() -  2.0 * (get_mla() 
		+  get_mra())) + 2.0 * (- 2.0 * kml * Lt * square(mm) - 2.0 * kt * Ls * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * square(mm) 
		+ 2.0 * kt * Ls * cos(get_bla() +  get_mra()) * square(mm) + 2.0 * kml * Lt * cos(2.0 * (get_mla() +  get_mra())) * square(mm) + 2.0 
		* kml * Ls * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) + 2.0 * kt * Ls * cos(get_bla() +  get_mla() +  2.0 * get_mra()) 
		* square(mm) - 2.0 * kml * Lt * mb * mm -  4.0 * kml * Lt * mt * mm - 2.0 * kt * Ls * mb * cos(get_bla() -  2.0 * get_mla() -  get_mra()) 
		* mm - 2.0 * kt * Ls * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm + 2.0 * kt * Ls * mb * cos(get_bla() +  get_mra()) 
		* mm + 2.0 * kt * Ls * mt * cos(get_bla() +  get_mra()) * mm + 2.0 * kml * Lt * mb * cos(2.0 * (get_mla() +  get_mra())) * mm + 2.0 
		* kml * Ls * mb * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * mm + 2.0 * kt * Ls * mb * cos(get_bla() +  get_mla() +  2.0 * get_mra()) 
		* mm - 2.0 * kml * Lt * mb * mt - 2.0 * (kml +  kt) * Ls * (mb * (mm +  mt) +  mm * (mm +  2.0 * mt)) * cos(get_bla() -  get_mla()) 
		+  (kml +  kt) * Ls * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) + 2.0 * kml * Lt * mb * mt * cos(2.0 
		* (get_bra() -  get_mra())) - kt * Ls * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mra()) - kt * Ls * mb * mt * cos(get_bla() 
		-  2.0 * get_mla() -  get_mra()) + kt * Ls * mb * mt * cos(get_bla() +  get_mra()) + kt * Ls * mb * mt * cos(get_bla() -  2.0 * get_bra() 
		-  2.0 * get_mla() +  get_mra()) + kml * Ls * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) + kt * Ls 
		* mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra())) * get_mla() -  2.0 * kmr * Lt * mb * mt * cos(get_bla() 
		+  get_bra()) * get_mra() + 2.0 * kmr * Lt * mb * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) * get_mra() - 4.0 * kt * Ls * 
		mb * mt * cos(get_bla() -  get_mla()) * get_mra() + 2.0 * kt * Ls * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 
		* get_mra()) * get_mra() -  2.0 * kmr * Ls * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mra()) * get_mra() - 2.0 * kt * Ls * 
		mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mra()) * get_mra() - 2.0 * kmr * Ls * mb * mt * cos(get_bla() -  2.0 * get_mla() 
		-  get_mra()) * get_mra() - 2.0 * kt * Ls * mb * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mra() + 2.0 * kmr * Ls * 
		mb * mt * cos(get_bla() +  get_mra()) * get_mra() + 2.0 * kt * Ls * mb * mt * cos(get_bla() +  get_mra()) * get_mra() + 2.0 * kmr * 
		Ls * mb * mt * cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() +  2.0 * kt * Ls * mb * mt * cos(get_bla() 
		-  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() + 2.0 * kmr * Lt * mb * mt * cos(get_bla() -  get_bra() +  2.0 * get_mra()) 
		* get_mra() + 2.0 * kt * Ls * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) * get_mra() -  2.0 * kmr 
		* Lt * mb * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * get_mra() - 4.0 * kmr * Lt * mb * mt * get_bra() * 
		sin(get_bla() -  get_mla()) * sin(get_bra() +  get_mla()) - 4.0 * kmr * Lt * mb * mt * get_bra() * sin(get_bla() -  get_mla()) * sin(get_bra() 
		-  get_mla() -  2.0 * get_mra()) - 4.0 * kmr * Ls * mb * mt * get_bra() * sin(get_bla() -  get_mla()) * sin(2.0 * get_bra() +  get_mla() 
		-  get_mra()) + 4.0 * kmr * Ls * mb * mt * get_bra() * sin(get_bla() -  get_mla()) * sin(get_mla() +  get_mra()) - 4.0 * kml * Lt * 
		mb * mt * M_PI)/(square(Ls) * Lt * (- 4.0 * cos(2.0 * (get_mla() +  get_mra())) * cube(mm) +  4.0 * cube(mm) +  4.0 * mb * square(mm) + 
		12.0 * mt * square(mm) -  4.0 * mb * cos(2.0 * (get_mla() +  get_mra())) * square(mm) + 4.0 * square(mt) * mm +  8.0 * mb * mt * mm 
		+  2.0 * mb * square(mt) - 2.0 * mt * (2.0 * mm * (mm +  mt) +  mb * (2.0 * mm +  mt)) * cos(2.0 * (get_bla() -  get_mla())) - 2.0 
		* mb * mt * (2.0 * mm +  mt) * cos(2.0 * (get_bra() -  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() +  get_bra() -  get_mla() 
		-  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() -  get_bra() -  get_mla() +  get_mra()))));
		
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_dd_mla() {
	if (state == LEFTBOUND) return 
		- (- 4.0 * square(Ls) * Lt * sin(get_bla() -  get_mla()) * square(get_d_bla()) * cube(mm) + 4.0 * square(Ls) * Lt * sin(get_bla() + get_mla()
		+ 2.0 * get_mra()) * square(get_d_bla()) * cube(mm) + 4.0 * Ls * square(Lt) * sin(2.0 * (get_mla() + get_mra())) * square(get_d_mla()) 
		* cube(mm) + 8.0 * Ls * square(Lt) * sin(get_mla() +  get_mra()) * square(get_d_mra()) * cube(mm) - 4.0 * square(Ls) * Lt * mb * sin(get_bla()
		-  get_mla()) * square(get_d_bla()) * square(mm) - 16.0 * square(Ls) * Lt * mt * sin(get_bla() -  get_mla()) * square(get_d_bla()) 
		* square(mm) + 4.0 * square(Ls) * Lt * mb * sin(get_bla() +  get_mla() +  2.0 * get_mra()) * square(get_d_bla()) * square(mm) + 4.0 
		* square(Ls) * Lt * mb * sin(get_bra() +  get_mla()) * square(get_d_bra()) * square(mm) - 4.0 * square(Ls) * Lt * mb * sin(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * square(get_d_bra()) * square(mm) - 4.0 * Ls * square(Lt) * mt * sin(2.0 * (get_bla() -  get_mla())) 
		* square(get_d_mla()) * square(mm) + 4.0 * Ls * square(Lt) * mb * sin(2.0 * (get_mla() +  get_mra())) * square(get_d_mla()) * square(mm) 
		- 4.0 * Ls * square(Lt) * mt * sin(2.0 * get_bla() -  get_mla() +  get_mra()) * square(get_d_mra()) * square(mm) + 8.0 * Ls * square(Lt) 
		* mb * sin(get_mla() +  get_mra()) * square(get_d_mra()) * square(mm) + 4.0 * Ls * square(Lt) * mt * sin(get_mla() +  get_mra()) * 
		square(get_d_mra()) * square(mm) - 12.0 * kml * Ls * mA * square(mm) +  12.0 * kt * Ls * tA * square(mm) + 4.0 * bA * kbl * Lt * cos(get_bla() 
		-  get_mla()) * square(mm) - 4.0 * kml * Lt * mA * cos(get_bla() -  get_mla()) * square(mm) + 4.0 * kml * Lt * M_PI * cos(get_bla() 
		-  get_mla()) * square(mm) + 4.0 * kmr * Lt * mA * cos(get_bra() +  get_mla()) * square(mm) - 4.0 * kmr * Lt * M_PI * cos(get_bra() 
		+  get_mla()) * square(mm) - 4.0 * kmr * Lt * mA * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * square(mm) + 4.0 * kmr * Lt * M_PI 
		* cos(get_bra() -  get_mla() -  2.0 * get_mra()) * square(mm) + 4.0 * kml * Ls * mA * cos(2.0 * (get_bla() +  get_mra())) * square(mm) 
		- 4.0 * kt * Ls * tA * cos(2.0 * (get_bla() +  get_mra())) * square(mm) - 4.0 * kml * Ls * M_PI * cos(2.0 * (get_bla() +  get_mra())) 
		* square(mm) + 4.0 * kt * Ls * M_PI * cos(2.0 * (get_bla() +  get_mra())) * square(mm) + 4.0 * kmr * Ls * mA * cos(2.0 * get_bla() 
		-  get_mla() +  get_mra()) * square(mm) - 4.0 * kt * Ls * tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * square(mm) - 4.0 * 
		kmr * Ls * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * square(mm) + 4.0 * kt * Ls * M_PI * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * square(mm) - 12.0 * kmr * Ls * mA * cos(get_mla() +  get_mra()) * square(mm) + 12.0 * kt * Ls * tA * cos(get_mla() 
		+  get_mra()) * square(mm) + 12.0 * kmr * Ls * M_PI * cos(get_mla() +  get_mra()) * square(mm) - 12.0 * kt * Ls * M_PI * cos(get_mla() 
		+  get_mra()) * square(mm) - 4.0 * bA * kbl * Lt * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) + 4.0 * kml * Lt * mA 
		* cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) - 4.0 * kml * Lt * M_PI * cos(get_bla() +  get_mla() +  2.0 * get_mra()) 
		* square(mm) + 12.0 * kml * Ls * get_mla() * square(mm) +  12.0 * kt * Ls * get_mla() * square(mm) + 4.0 * kml * Lt * cos(get_bla() 
		-  get_mla()) * get_mla() * square(mm) - 4.0 * kml * Ls * cos(2.0 * (get_bla() +  get_mra())) * get_mla() * square(mm) - 4.0 * kt * 
		Ls * cos(2.0 * (get_bla() +  get_mra())) * get_mla() * square(mm) - 4.0 * kt * Ls * cos(2.0 * get_bla() -  get_mla() +  get_mra()) 
		* get_mla() * square(mm) + 12.0 * kt * Ls * cos(get_mla() +  get_mra()) * get_mla() * square(mm) - 4.0 * kml * Lt * cos(get_bla() + 
		 get_mla() +  2.0 * get_mra()) * get_mla() * square(mm) + 12.0 * kt * Ls * get_mra() * square(mm) - 4.0 * kmr * Lt * cos(get_bra() 
		+  get_mla()) * get_mra() * square(mm) + 4.0 * kmr * Lt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mra() * square(mm) - 
		4.0 * kt * Ls * cos(2.0 * (get_bla() +  get_mra())) * get_mra() * square(mm) - 4.0 * kmr * Ls * cos(2.0 * get_bla() -  get_mla() + 
		 get_mra()) * get_mra() * square(mm) - 4.0 * kt * Ls * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() * square(mm) + 12.0 
		* kmr * Ls * cos(get_mla() +  get_mra()) * get_mra() * square(mm) + 12.0 * kt * Ls * cos(get_mla() +  get_mra()) * get_mra() * square(mm) 
		+ 12.0 * kml * Ls * M_PI * square(mm) -  12.0 * kt * Ls * M_PI * square(mm) - 8.0 * square(Ls) * Lt * square(mt) * sin(get_bla() - 
		 get_mla()) * square(get_d_bla()) * mm - 12.0 * square(Ls) * Lt * mb * mt * sin(get_bla() -  get_mla()) * square(get_d_bla()) * mm 
		+ 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) * square(get_d_bla()) * mm + 
		2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) * square(get_d_bla()) * mm - 2.0 
		* square(Ls) * Lt * mb * mt * sin(2.0 * get_bla() +  get_bra() -  get_mla()) * square(get_d_bra()) * mm + 2.0 * square(Ls) * Lt * mb 
		* mt * sin(get_bra() +  get_mla()) * square(get_d_bra()) * mm - 2.0 * square(Ls) * Lt * mb * mt * sin(get_bra() -  get_mla() -  2.0 
		* get_mra()) * square(get_d_bra()) * mm - 2.0 * square(Ls) * Lt * mb * mt * sin(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * 
		get_mra()) * square(get_d_bra()) * mm - 4.0 * Ls * square(Lt) * square(mt) * sin(2.0 * (get_bla() -  get_mla())) * square(get_d_mla()) 
		* mm - 4.0 * Ls * square(Lt) * mb * mt * sin(2.0 * (get_bla() -  get_mla())) * square(get_d_mla()) * mm - 4.0 * Ls * square(Lt) * mb 
		* mt * sin(2.0 * get_bla() -  get_mla() +  get_mra()) * square(get_d_mra()) * mm + 4.0 * Ls * square(Lt) * mb * mt * sin(get_mla() 
		+  get_mra()) * square(get_d_mra()) * mm -  8.0 * kml * Ls * mA * mb * mm -  8.0 * kml * Ls * mA * mt * mm + 8.0 * kt * Ls * mb * tA 
		* mm +  8.0 * kt * Ls * mt * tA * mm + 4.0 * bA * kbl * Lt * mb * cos(get_bla() -  get_mla()) * mm - 4.0 * kml * Lt * mA * mb * cos(get_bla() 
		-  get_mla()) * mm + 8.0 * bA * kbl * Lt * mt * cos(get_bla() -  get_mla()) * mm - 8.0 * kml * Lt * mA * mt * cos(get_bla() -  get_mla()) 
		* mm + 4.0 * kml * Lt * mb * M_PI * cos(get_bla() -  get_mla()) * mm + 8.0 * kml * Lt * mt * M_PI * cos(get_bla() -  get_mla()) * mm 
		- 2.0 * kmr * Lt * mA * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * mm + 2.0 * kmr * Lt * mt * M_PI * cos(2.0 * get_bla() 
		+  get_bra() -  get_mla()) * mm + 4.0 * kmr * Lt * mA * mb * cos(get_bra() +  get_mla()) * mm + 2.0 * kmr * Lt * mA * mt * cos(get_bra() 
		+  get_mla()) * mm - 4.0 * kmr * Lt * mb * M_PI * cos(get_bra() +  get_mla()) * mm - 2.0 * kmr * Lt * mt * M_PI * cos(get_bra() +  
		get_mla()) * mm - 4.0 * kmr * Lt * mA * mb * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm - 2.0 * kmr * Lt * mA * mt * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * mm + 4.0 * kmr * Lt * mb * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm + 2.0 * kmr 
		* Lt * mt * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm + 4.0 * kml * Ls * mA * mb * cos(2.0 * (get_bra() -  get_mra())) 
		* mm - 4.0 * kt * Ls * mb * tA * cos(2.0 * (get_bra() -  get_mra())) * mm - 4.0 * kml * Ls * mb * M_PI * cos(2.0 * (get_bra() -  get_mra())) 
		* mm + 4.0 * kt * Ls * mb * M_PI * cos(2.0 * (get_bra() -  get_mra())) * mm + 4.0 * kmr * Ls * mA * mb * cos(2.0 * get_bra() +  get_mla() 
		-  get_mra()) * mm - 4.0 * kt * Ls * mb * tA * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm - 4.0 * kmr * Ls * mb * M_PI * cos(2.0 
		* get_bra() +  get_mla() -  get_mra()) * mm + 4.0 * kt * Ls * mb * M_PI * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm + 4.0 
		* kml * Ls * mA * mb * cos(2.0 * (get_bla() +  get_mra())) * mm - 4.0 * kt * Ls * mb * tA * cos(2.0 * (get_bla() +  get_mra())) * mm 
		- 4.0 * kml * Ls * mb * M_PI * cos(2.0 * (get_bla() +  get_mra())) * mm + 4.0 * kt * Ls * mb * M_PI * cos(2.0 * (get_bla() +  get_mra())) 
		* mm + 4.0 * kmr * Ls * mA * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm + 4.0 * kmr * Ls * mA * mt * cos(2.0 * get_bla() 
		-  get_mla() +  get_mra()) * mm - 4.0 * kt * Ls * mb * tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kt * Ls * mt 
		* tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kmr * Ls * mb * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) 
		* mm + 4.0 * kt * Ls * mb * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kmr * Ls * mt * M_PI * cos(2.0 * get_bla() 
		-  get_mla() +  get_mra()) * mm + 4.0 * kt * Ls * mt * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 8.0 * kmr * Ls 
		* mA * mb * cos(get_mla() +  get_mra()) * mm - 4.0 * kmr * Ls * mA * mt * cos(get_mla() +  get_mra()) * mm + 8.0 * kt * Ls * mb * tA 
		* cos(get_mla() +  get_mra()) * mm + 4.0 * kt * Ls * mt * tA * cos(get_mla() +  get_mra()) * mm + 8.0 * kmr * Ls * mb * M_PI * cos(get_mla() 
		+  get_mra()) * mm - 8.0 * kt * Ls * mb * M_PI * cos(get_mla() +  get_mra()) * mm + 4.0 * kmr * Ls * mt * M_PI * cos(get_mla() +  get_mra()) 
		* mm - 4.0 * kt * Ls * mt * M_PI * cos(get_mla() +  get_mra()) * mm + 2.0 * kmr * Lt * mA * mt * cos(2.0 * get_bla() -  get_bra() - 
		 get_mla() +  2.0 * get_mra()) * mm - 2.0 * kmr * Lt * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) 
		* mm - 4.0 * bA * kbl * Lt * mb * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * mm + 4.0 * kml * Lt * mA * mb * cos(get_bla() + 
		 get_mla() +  2.0 * get_mra()) * mm - 4.0 * kml * Lt * mb * M_PI * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * mm + 8.0 * kml 
		* Ls * mb * get_mla() * mm +  8.0 * kt * Ls * mb * get_mla() * mm + 8.0 * kml * Ls * mt * get_mla() * mm +  8.0 * kt * Ls * mt * get_mla() 
		* mm + 4.0 * kml * Lt * mb * cos(get_bla() -  get_mla()) * get_mla() * mm + 8.0 * kml * Lt * mt * cos(get_bla() -  get_mla()) * get_mla() 
		* mm - 4.0 * kml * Ls * mb * cos(2.0 * (get_bra() -  get_mra())) * get_mla() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bra() -  get_mra())) 
		* get_mla() * mm - 4.0 * kt * Ls * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mla() * mm - 4.0 * kml * Ls * mb * cos(2.0 
		* (get_bla() +  get_mra())) * get_mla() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bla() +  get_mra())) * get_mla() * mm - 4.0 * kt 
		* Ls * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() * mm - 4.0 * kt * Ls * mt * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * get_mla() * mm + 8.0 * kt * Ls * mb * cos(get_mla() +  get_mra()) * get_mla() * mm + 4.0 * kt * Ls * mt * cos(get_mla() 
		+  get_mra()) * get_mla() * mm - 4.0 * kml * Lt * mb * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * get_mla() * mm + 8.0 * kt * 
		Ls * mb * get_mra() * mm +  8.0 * kt * Ls * mt * get_mra() * mm + 2.0 * kmr * Lt * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) 
		* get_mra() * mm - 4.0 * kmr * Lt * mb * cos(get_bra() +  get_mla()) * get_mra() * mm - 2.0 * kmr * Lt * mt * cos(get_bra() +  get_mla()) 
		* get_mra() * mm + 4.0 * kmr * Lt * mb * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mra() * mm + 2.0 * kmr * Lt * mt * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bra() -  get_mra())) * get_mra() * mm - 4.0 
		* kmr * Ls * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * get_bra() +  get_mla() 
		-  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bla() +  get_mra())) * get_mra() * mm - 4.0 * kmr * Ls * mb * 
		cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) 
		* get_mra() * mm - 4.0 * kmr * Ls * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mt * cos(2.0 
		* get_bla() -  get_mla() +  get_mra()) * get_mra() * mm + 8.0 * kmr * Ls * mb * cos(get_mla() +  get_mra()) * get_mra() * mm + 8.0 
		* kt * Ls * mb * cos(get_mla() +  get_mra()) * get_mra() * mm + 4.0 * kmr * Ls * mt * cos(get_mla() +  get_mra()) * get_mra() * mm 
		+ 4.0 * kt * Ls * mt * cos(get_mla() +  get_mra()) * get_mra() * mm - 2.0 * kmr * Lt * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() 
		+  2.0 * get_mra()) * get_mra() * mm +  8.0 * kml * Ls * mb * M_PI * mm -  8.0 * kt * Ls * mb * M_PI * mm + 8.0 * kml * Ls * mt * M_PI 
		* mm -  8.0 * kt * Ls * mt * M_PI * mm - 4.0 * square(Ls) * Lt * mb * square(mt) * sin(get_bla() -  get_mla()) * square(get_d_bla()) 
		+  2.0 * square(Ls) * Lt * mb * square(mt) * sin(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) * square(get_d_bla()) 
		+ 2.0 * square(Ls) * Lt * mb * square(mt) * sin(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) * square(get_d_bla()) 
		-  2.0 * Ls * square(Lt) * mb * square(mt) * sin(2.0 * (get_bla() -  get_mla())) * square(get_d_mla()) + Ls * square(Lt) * mb * square(mt) 
		* sin(2.0 * (get_bla() +  get_bra() -  get_mla() -  get_mra())) * square(get_d_mla()) +  Ls * square(Lt) * mb * square(mt) * sin(2.0 
		* (get_bla() -  get_bra() -  get_mla() +  get_mra())) * square(get_d_mla()) -  4.0 * kml * Ls * mA * mb * mt +  4.0 * kt * Ls * mb 
		* mt * tA + 4.0 * bA * kbl * Lt * mb * mt * cos(get_bla() -  get_mla()) - 4.0 * kml * Lt * mA * mb * mt * cos(get_bla() -  get_mla()) 
		+ 4.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() -  get_mla()) - 2.0 * kmr * Lt * mA * mb * mt * cos(2.0 * get_bla() +  get_bra() 
		-  get_mla()) + 2.0 * kmr * Lt * mb * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  get_mla()) + 2.0 * kmr * Lt * mA * mb * mt * 
		cos(get_bra() +  get_mla()) - 2.0 * kmr * Lt * mb * mt * M_PI * cos(get_bra() +  get_mla()) - 2.0 * kmr * Lt * mA * mb * mt * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) + 2.0 * kmr * Lt * mb * mt * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) - 2.0 * bA * kbl 
		* Lt * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) + 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() 
		+  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) - 2.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() +  2.0 * get_bra() -  get_mla() 
		-  2.0 * get_mra()) + 4.0 * kml * Ls * mA * mb * mt * cos(2.0 * (get_bra() -  get_mra())) - 4.0 * kt * Ls * mb * mt * tA * cos(2.0 
		* (get_bra() -  get_mra())) - 4.0 * kml * Ls * mb * mt * M_PI * cos(2.0 * (get_bra() -  get_mra())) + 4.0 * kt * Ls * mb * mt * M_PI 
		* cos(2.0 * (get_bra() -  get_mra())) - 2.0 * kmr * Ls * mA * mb * mt * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) 
		+ 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) + 2.0 * kmr * Ls * mb * mt * M_PI 
		* cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) - 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bla() +  2.0 * 
		get_bra() -  get_mla() -  get_mra()) + 2.0 * kmr * Ls * mA * mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) - 2.0 * kt * 
		Ls * mb * mt * tA * cos(2.0 * get_bra() +  get_mla() -  get_mra()) - 2.0 * kmr * Ls * mb * mt * M_PI * cos(2.0 * get_bra() +  get_mla() 
		-  get_mra()) + 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bra() +  get_mla() -  get_mra()) + 2.0 * kmr * Ls * mA * mb * mt * cos(2.0 
		* get_bla() -  get_mla() +  get_mra()) - 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) - 2.0 * kmr 
		* Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) + 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) - 2.0 * kmr * Ls * mA * mb * mt * cos(get_mla() +  get_mra()) + 2.0 * kt * Ls * mb * mt * tA * cos(get_mla() +  get_mra()) 
		+ 2.0 * kmr * Ls * mb * mt * M_PI * cos(get_mla() +  get_mra()) - 2.0 * kt * Ls * mb * mt * M_PI * cos(get_mla() +  get_mra()) - 2.0 
		* bA * kbl * Lt * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) + 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() 
		-  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) - 2.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() -  2.0 * get_bra() -  get_mla() 
		+  2.0 * get_mra()) + 2.0 * kmr * Lt * mA * mb * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) - 2.0 * kmr 
		* Lt * mb * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) + 2.0 * kmr * get_bra() * (- 2.0 * Lt * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * square(mm) + 2.0 * Ls * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * square(mm) - 6.0 * Ls * 
		cos(get_mla() +  get_mra()) * square(mm) - 2.0 * Lt * mb * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm - Lt * mt * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * mm + 2.0 * Ls * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm + 2.0 * Ls * mb * cos(2.0 
		* get_bla() -  get_mla() +  get_mra()) * mm + 2.0 * Ls * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * Ls * mb * 
		cos(get_mla() +  get_mra()) * mm - 2.0 * Ls * mt * cos(get_mla() +  get_mra()) * mm + Lt * mt * cos(2.0 * get_bla() -  get_bra() - 
		 get_mla() +  2.0 * get_mra()) * mm - Lt * (mb +  mm) * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) + Lt * (mb +  mm) * (2.0 
		* mm +  mt) * cos(get_bra() +  get_mla()) - Lt * mb * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) - Ls * mb * mt * cos(2.0 
		* get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) + Ls * mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) + Ls * mb 
		* mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) - Ls * mb * mt * cos(get_mla() +  get_mra()) + Lt * mb * mt * cos(2.0 * get_bla() 
		-  get_bra() -  get_mla() +  2.0 * get_mra())) + 2.0 * get_bla() * (- 6.0 * kml * Ls * square(mm) + 2.0 * kml * Ls * cos(2.0 * (get_bla() 
		+  get_mra())) * square(mm) + 2.0 * kbl * Lt * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * square(mm) + 2.0 * kml * Lt * cos(get_bla() 
		+  get_mla() +  2.0 * get_mra()) * square(mm) - 4.0 * kml * Ls * mb * mm -  4.0 * kml * Ls * mt * mm + 2.0 * kml * Ls * mb * cos(2.0 
		* (get_bra() -  get_mra())) * mm + 2.0 * kml * Ls * mb * cos(2.0 * (get_bla() +  get_mra())) * mm + 2.0 * kbl * Lt * mb * cos(get_bla() 
		+  get_mla() +  2.0 * get_mra()) * mm + 2.0 * kml * Lt * mb * cos(get_bla() +  get_mla() +  2.0 * get_mra()) * mm - 2.0 * kml * Ls 
		* mb * mt - 2.0 * (kbl +  kml) * Lt * (mb * (mm +  mt) +  mm * (mm +  2.0 * mt)) * cos(get_bla() -  get_mla()) +  (kbl +  kml) * Lt 
		* mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mla() -  2.0 * get_mra()) + 2.0 * kml * Ls * mb * mt * cos(2.0 * (get_bra() -  
		get_mra())) + kbl * Lt * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) + kml * Lt * mb * mt * cos(get_bla() 
		-  2.0 * get_bra() -  get_mla() +  2.0 * get_mra())) + 4.0 * kml * Ls * mb * mt * get_mla() +  4.0 * kt * Ls * mb * mt * get_mla() 
		+ 4.0 * kml * Lt * mb * mt * cos(get_bla() -  get_mla()) * get_mla() - 2.0 * kml * Lt * mb * mt * cos(get_bla() +  2.0 * get_bra() 
		-  get_mla() -  2.0 * get_mra()) * get_mla() -  4.0 * kml * Ls * mb * mt * cos(2.0 * (get_bra() -  get_mra())) * get_mla() - 4.0 * 
		kt * Ls * mb * mt * cos(2.0 * (get_bra() -  get_mra())) * get_mla() + 2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() +  2.0 * get_bra() 
		-  get_mla() -  get_mra()) * get_mla() -  2.0 * kt * Ls * mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mla() - 2.0 
		* kt * Ls * mb * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() + 2.0 * kt * Ls * mb * mt * cos(get_mla() +  get_mra()) 
		* get_mla() - 2.0 * kml * Lt * mb * mt * cos(get_bla() -  2.0 * get_bra() -  get_mla() +  2.0 * get_mra()) * get_mla() +  4.0 * kt 
		* Ls * mb * mt * get_mra() + 2.0 * kmr * Lt * mb * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * get_mra() - 2.0 * kmr * Lt 
		* mb * mt * cos(get_bra() +  get_mla()) * get_mra() + 2.0 * kmr * Lt * mb * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mra() 
		- 4.0 * kt * Ls * mb * mt * cos(2.0 * (get_bra() -  get_mra())) * get_mra() + 2.0 * kmr * Ls * mb * mt * cos(2.0 * get_bla() +  2.0 
		* get_bra() -  get_mla() -  get_mra()) * get_mra() +  2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() 
		-  get_mra()) * get_mra() - 2.0 * kmr * Ls * mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mra() - 2.0 * kt * Ls * 
		mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mra() - 2.0 * kmr * Ls * mb * mt * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * get_mra() - 2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() + 2.0 * kmr * Ls * 
		mb * mt * cos(get_mla() +  get_mra()) * get_mra() + 2.0 * kt * Ls * mb * mt * cos(get_mla() +  get_mra()) * get_mra() - 2.0 * kmr * 
		Lt * mb * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * get_mra() +  4.0 * kml * Ls * mb * mt * M_PI - 4.0 
		* kt * Ls * mb * mt * M_PI)/(Ls * square(Lt) * (- 4.0 * cos(2.0 * (get_mla() +  get_mra())) * cube(mm) +  4.0 * cube(mm) +  4.0 * mb * 
		square(mm) + 12.0 * mt * square(mm) -  4.0 * mb * cos(2.0 * (get_mla() +  get_mra())) * square(mm) + 4.0 * square(mt) * mm +  8.0 * 
		mb * mt * mm +  2.0 * mb * square(mt) - 2.0 * mt * (2.0 * mm * (mm +  mt) +  mb * (2.0 * mm +  mt)) * cos(2.0 * (get_bla() -  get_mla())) 
		- 2.0 * mb * mt * (2.0 * mm +  mt) * cos(2.0 * (get_bra() -  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() +  get_bra() -  get_mla() 
		-  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() -  get_bra() -  get_mla() +  get_mra()))));
		
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_dd_mra() {
	if (state == LEFTBOUND) return 
		- (- 4.0 * square(Ls) * Lt * sin(get_bla() -  2.0 * get_mla() -  get_mra()) * square(get_d_bla()) * cube(mm) + 4.0 * square(Ls) * Lt 
		* sin(get_bla() +  get_mra()) * square(get_d_bla()) * cube(mm) + 8.0 * Ls * square(Lt) * sin(get_mla() +  get_mra()) * square(get_d_mla()) 
		* cube(mm) + 4.0 * Ls * square(Lt) * sin(2.0 * (get_mla() +  get_mra())) * square(get_d_mra()) * cube(mm) - 4.0 * square(Ls) * Lt * mb 
		* sin(get_bla() -  2.0 * get_mla() -  get_mra()) * square(get_d_bla()) * square(mm) - 4.0 * square(Ls) * Lt * mt * sin(get_bla() - 
		2.0 * get_mla() - get_mra()) * square(get_d_bla()) * square(mm) + 4.0 * square(Ls) * Lt * mb * sin(get_bla() +  get_mra()) * square(get_d_bla()) 
		* square(mm) + 4.0 * square(Ls) * Lt * mt * sin(get_bla() +  get_mra()) * square(get_d_bla()) * square(mm) - 4.0 * square(Ls) * Lt 
		* mb * sin(get_bra() -  get_mra()) * square(get_d_bra()) * square(mm) + 4.0 * square(Ls) * Lt * mb * sin(get_bra() +  2.0 * get_mla() 
		+  get_mra()) * square(get_d_bra()) * square(mm) + 8.0 * Ls * square(Lt) * mb * sin(get_mla() +  get_mra()) * square(get_d_mla()) * 
		square(mm) + 8.0 * Ls * square(Lt) * mt * sin(get_mla() +  get_mra()) * square(get_d_mla()) * square(mm) + 4.0 * Ls * square(Lt) * 
		mb * sin(2.0 * (get_mla() +  get_mra())) * square(get_d_mra()) * square(mm) -  12.0 * kmr * Ls * mA * square(mm) +  12.0 * kt * Ls 
		* tA * square(mm) + 4.0 * kmr * Ls * mA * cos(2.0 * (get_bla() -  get_mla())) * square(mm) - 4.0 * kt * Ls * tA * cos(2.0 * (get_bla() 
		-  get_mla())) * square(mm) - 4.0 * kmr * Ls * M_PI * cos(2.0 * (get_bla() -  get_mla())) * square(mm) + 4.0 * kt * Ls * M_PI * cos(2.0 
		* (get_bla() -  get_mla())) * square(mm) - 4.0 * kmr * Lt * mA * cos(get_bra() -  get_mra()) * square(mm) + 4.0 * kmr * Lt * M_PI * 
		cos(get_bra() -  get_mra()) * square(mm) + 4.0 * bA * kbl * Lt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * square(mm) - 4.0 
		* kml * Lt * mA * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * square(mm) + 4.0 * kml * Lt * M_PI * cos(get_bla() -  2.0 * get_mla() 
		-  get_mra()) * square(mm) - 4.0 * bA * kbl * Lt * cos(get_bla() +  get_mra()) * square(mm) + 4.0 * kml * Lt * mA * cos(get_bla() + 
		 get_mra()) * square(mm) - 4.0 * kml * Lt * M_PI * cos(get_bla() +  get_mra()) * square(mm) + 4.0 * kml * Ls * mA * cos(2.0 * get_bla() 
		-  get_mla() +  get_mra()) * square(mm) - 4.0 * kt * Ls * tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * square(mm) - 4.0 * 
		kml * Ls * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * square(mm) + 4.0 * kt * Ls * M_PI * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * square(mm) - 12.0 * kml * Ls * mA * cos(get_mla() +  get_mra()) * square(mm) + 12.0 * kt * Ls * tA * cos(get_mla() 
		+  get_mra()) * square(mm) + 12.0 * kml * Ls * M_PI * cos(get_mla() +  get_mra()) * square(mm) - 12.0 * kt * Ls * M_PI * cos(get_mla() 
		+  get_mra()) * square(mm) + 4.0 * kmr * Lt * mA * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * square(mm) - 4.0 * kmr * Lt * M_PI 
		* cos(get_bra() +  2.0 * get_mla() +  get_mra()) * square(mm) + 12.0 * kt * Ls * get_mla() * square(mm) - 4.0 * kt * Ls * cos(2.0 * 
		(get_bla() -  get_mla())) * get_mla() * square(mm) + 4.0 * kml * Lt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mla() * 
		square(mm) - 4.0 * kml * Lt * cos(get_bla() +  get_mra()) * get_mla() * square(mm) - 4.0 * kml * Ls * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * get_mla() * square(mm) - 4.0 * kt * Ls * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() * square(mm) + 
		12.0 * kml * Ls * cos(get_mla() +  get_mra()) * get_mla() * square(mm) + 12.0 * kt * Ls * cos(get_mla() +  get_mra()) * get_mla() * 
		square(mm) + 12.0 * kmr * Ls * get_mra() * square(mm) +  12.0 * kt * Ls * get_mra() * square(mm) - 4.0 * kmr * Ls * cos(2.0 * (get_bla() 
		-  get_mla())) * get_mra() * square(mm) - 4.0 * kt * Ls * cos(2.0 * (get_bla() -  get_mla())) * get_mra() * square(mm) + 4.0 * kmr 
		* Lt * cos(get_bra() -  get_mra()) * get_mra() * square(mm) - 4.0 * kt * Ls * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() 
		* square(mm) + 12.0 * kt * Ls * cos(get_mla() +  get_mra()) * get_mra() * square(mm) - 4.0 * kmr * Lt * cos(get_bra() +  2.0 * get_mla() 
		+  get_mra()) * get_mra() * square(mm) + 12.0 * kmr * Ls * M_PI * square(mm) -  12.0 * kt * Ls * M_PI * square(mm) - 2.0 * square(Ls) 
		* Lt * mb * mt * sin(get_bla() +  2.0 * get_bra() -  get_mra()) * square(get_d_bla()) * mm - 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() 
		-  2.0 * get_mla() -  get_mra()) * square(get_d_bla()) * mm + 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() +  get_mra()) * square(get_d_bla()) 
		* mm + 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) * square(get_d_bla()) * 
		mm - 12.0 * square(Ls) * Lt * mb * mt * sin(get_bra() -  get_mra()) * square(get_d_bra()) * mm + 2.0 * square(Ls) * Lt * mb * mt * 
		sin(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * square(get_d_bra()) * mm - 2.0 * square(Ls) * Lt * mb * mt * sin(2.0 
		* get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * square(get_d_bra()) * mm - 4.0 * Ls * square(Lt) * mb * mt * sin(2.0 * 
		get_bra() +  get_mla() -  get_mra()) * square(get_d_mla()) * mm + 4.0 * Ls * square(Lt) * mb * mt * sin(get_mla() +  get_mra()) * square(get_d_mla()) 
		* mm - 4.0 * Ls * square(Lt) * mb * mt * sin(2.0 * (get_bra() -  get_mra())) * square(get_d_mra()) * mm -  8.0 * kmr * Ls * mA * mb 
		* mm -  16.0 * kmr * Ls * mA * mt * mm + 8.0 * kt * Ls * mb * tA * mm +  16.0 * kt * Ls * mt * tA * mm + 4.0 * kmr * Ls * mA * mb * 
		cos(2.0 * (get_bla() -  get_mla())) * mm + 8.0 * kmr * Ls * mA * mt * cos(2.0 * (get_bla() -  get_mla())) * mm - 4.0 * kt * Ls * mb 
		* tA * cos(2.0 * (get_bla() -  get_mla())) * mm - 8.0 * kt * Ls * mt * tA * cos(2.0 * (get_bla() -  get_mla())) * mm - 4.0 * kmr * 
		Ls * mb * M_PI * cos(2.0 * (get_bla() -  get_mla())) * mm + 4.0 * kt * Ls * mb * M_PI * cos(2.0 * (get_bla() -  get_mla())) * mm - 
		8.0 * kmr * Ls * mt * M_PI * cos(2.0 * (get_bla() -  get_mla())) * mm + 8.0 * kt * Ls * mt * M_PI * cos(2.0 * (get_bla() -  get_mla())) 
		* mm + 4.0 * kmr * Ls * mA * mb * cos(2.0 * (get_bra() +  get_mla())) * mm - 4.0 * kt * Ls * mb * tA * cos(2.0 * (get_bra() +  get_mla())) 
		* mm - 4.0 * kmr * Ls * mb * M_PI * cos(2.0 * (get_bra() +  get_mla())) * mm + 4.0 * kt * Ls * mb * M_PI * cos(2.0 * (get_bra() +  
		get_mla())) * mm - 4.0 * kmr * Lt * mA * mb * cos(get_bra() -  get_mra()) * mm - 12.0 * kmr * Lt * mA * mt * cos(get_bra() -  get_mra()) 
		* mm + 4.0 * kmr * Lt * mb * M_PI * cos(get_bra() -  get_mra()) * mm + 12.0 * kmr * Lt * mt * M_PI * cos(get_bra() -  get_mra()) * 
		mm + 4.0 * bA * kbl * Lt * mb * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm - 4.0 * kml * Lt * mA * mb * cos(get_bla() -  2.0 
		* get_mla() -  get_mra()) * mm + 4.0 * bA * kbl * Lt * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm - 4.0 * kml * Lt * 
		mA * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm + 4.0 * kml * Lt * mb * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) 
		* mm + 4.0 * kml * Lt * mt * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * mm + 2.0 * kmr * Lt * mA * mt * cos(2.0 * get_bla() 
		+  get_bra() -  2.0 * get_mla() -  get_mra()) * mm - 2.0 * kmr * Lt * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() 
		-  get_mra()) * mm + 4.0 * kml * Ls * mA * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm - 4.0 * kt * Ls * mb * tA * cos(2.0 
		* get_bra() +  get_mla() -  get_mra()) * mm - 4.0 * kml * Ls * mb * M_PI * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm + 4.0 
		* kt * Ls * mb * M_PI * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm - 4.0 * bA * kbl * Lt * mb * cos(get_bla() +  get_mra()) 
		* mm + 4.0 * kml * Lt * mA * mb * cos(get_bla() +  get_mra()) * mm - 4.0 * bA * kbl * Lt * mt * cos(get_bla() +  get_mra()) * mm + 
		4.0 * kml * Lt * mA * mt * cos(get_bla() +  get_mra()) * mm - 4.0 * kml * Lt * mb * M_PI * cos(get_bla() +  get_mra()) * mm - 4.0 * 
		kml * Lt * mt * M_PI * cos(get_bla() +  get_mra()) * mm + 2.0 * kmr * Lt * mA * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() 
		+  get_mra()) * mm - 2.0 * kmr * Lt * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * mm + 4.0 * kml 
		* Ls * mA * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm + 4.0 * kml * Ls * mA * mt * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * mm - 4.0 * kt * Ls * mb * tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kt * Ls * mt * tA * cos(2.0 
		* get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kml * Ls * mb * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm + 4.0 
		* kt * Ls * mb * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kml * Ls * mt * M_PI * cos(2.0 * get_bla() -  get_mla() 
		+  get_mra()) * mm + 4.0 * kt * Ls * mt * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 8.0 * kml * Ls * mA * mb * cos(get_mla() 
		+  get_mra()) * mm - 4.0 * kml * Ls * mA * mt * cos(get_mla() +  get_mra()) * mm + 8.0 * kt * Ls * mb * tA * cos(get_mla() +  get_mra()) 
		* mm + 4.0 * kt * Ls * mt * tA * cos(get_mla() +  get_mra()) * mm + 8.0 * kml * Ls * mb * M_PI * cos(get_mla() +  get_mra()) * mm - 
		8.0 * kt * Ls * mb * M_PI * cos(get_mla() +  get_mra()) * mm + 4.0 * kml * Ls * mt * M_PI * cos(get_mla() +  get_mra()) * mm - 4.0 
		* kt * Ls * mt * M_PI * cos(get_mla() +  get_mra()) * mm + 4.0 * kmr * Lt * mA * mb * cos(get_bra() +  2.0 * get_mla() +  get_mra()) 
		* mm - 4.0 * kmr * Lt * mb * M_PI * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * mm + 8.0 * kt * Ls * mb * get_mla() * mm +  16.0 
		* kt * Ls * mt * get_mla() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bla() -  get_mla())) * get_mla() * mm - 8.0 * kt * Ls * mt * 
		cos(2.0 * (get_bla() -  get_mla())) * get_mla() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bra() +  get_mla())) * get_mla() * mm + 
		4.0 * kml * Lt * mb * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mla() * mm + 4.0 * kml * Lt * mt * cos(get_bla() -  2.0 
		* get_mla() -  get_mra()) * get_mla() * mm - 4.0 * kml * Ls * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mla() * mm 
		- 4.0 * kt * Ls * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mla() * mm - 4.0 * kml * Lt * mb * cos(get_bla() +  get_mra()) 
		* get_mla() * mm - 4.0 * kml * Lt * mt * cos(get_bla() +  get_mra()) * get_mla() * mm - 4.0 * kml * Ls * mb * cos(2.0 * get_bla() - 
		 get_mla() +  get_mra()) * get_mla() * mm - 4.0 * kt * Ls * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() * mm - 
		4.0 * kml * Ls * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() * mm - 4.0 * kt * Ls * mt * cos(2.0 * get_bla() - 
		 get_mla() +  get_mra()) * get_mla() * mm + 8.0 * kml * Ls * mb * cos(get_mla() +  get_mra()) * get_mla() * mm + 8.0 * kt * Ls * mb 
		* cos(get_mla() +  get_mra()) * get_mla() * mm + 4.0 * kml * Ls * mt * cos(get_mla() +  get_mra()) * get_mla() * mm + 4.0 * kt * Ls 
		* mt * cos(get_mla() +  get_mra()) * get_mla() * mm + 8.0 * kmr * Ls * mb * get_mra() * mm +  8.0 * kt * Ls * mb * get_mra() * mm + 
		16.0 * kmr * Ls * mt * get_mra() * mm +  16.0 * kt * Ls * mt * get_mra() * mm - 4.0 * kmr * Ls * mb * cos(2.0 * (get_bla() -  get_mla())) 
		* get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bla() -  get_mla())) * get_mra() * mm - 8.0 * kmr * Ls * mt * cos(2.0 * (get_bla() 
		-  get_mla())) * get_mra() * mm - 8.0 * kt * Ls * mt * cos(2.0 * (get_bla() -  get_mla())) * get_mra() * mm - 4.0 * kmr * Ls * mb * 
		cos(2.0 * (get_bra() +  get_mla())) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * (get_bra() +  get_mla())) * get_mra() * mm + 
		4.0 * kmr * Lt * mb * cos(get_bra() -  get_mra()) * get_mra() * mm + 12.0 * kmr * Lt * mt * cos(get_bra() -  get_mra()) * get_mra() 
		* mm - 2.0 * kmr * Lt * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() * mm -  4.0 * kt * Ls * 
		mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mra() * mm -  2.0 * kmr * Lt * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 
		* get_mla() +  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() * mm - 
		4.0 * kt * Ls * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() * mm + 8.0 * kt * Ls * mb * cos(get_mla() +  get_mra()) 
		* get_mra() * mm + 4.0 * kt * Ls * mt * cos(get_mla() +  get_mra()) * get_mra() * mm - 4.0 * kmr * Lt * mb * cos(get_bra() +  2.0 * 
		get_mla() +  get_mra()) * get_mra() * mm + 8.0 * kmr * Ls * mb * M_PI * mm -  8.0 * kt * Ls * mb * M_PI * mm + 16.0 * kmr * Ls * mt 
		* M_PI * mm -  16.0 * kt * Ls * mt * M_PI * mm - 4.0 * kmr * Ls * mA * square(mt) - 4.0 * square(Ls) * Lt * mb * square(mt) * sin(get_bra() 
		-  get_mra()) * square(get_d_bra()) + 2.0 * square(Ls) * Lt * mb * square(mt) * sin(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() 
		-  get_mra()) * square(get_d_bra()) -  2.0 * square(Ls) * Lt * mb * square(mt) * sin(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() 
		+  get_mra()) * square(get_d_bra()) -  2.0 * Ls * square(Lt) * mb * square(mt) * sin(2.0 * (get_bra() -  get_mra())) * square(get_d_mra()) 
		+ Ls * square(Lt) * mb * square(mt) * sin(2.0 * (get_bla() +  get_bra() -  get_mla() -  get_mra())) * square(get_d_mra()) -  Ls * square(Lt) 
		* mb * square(mt) * sin(2.0 * (get_bla() -  get_bra() -  get_mla() +  get_mra())) * square(get_d_mra()) -  4.0 * kmr * Ls * mA * mb 
		* mt +  4.0 * kt * Ls * square(mt) * tA + 4.0 * kt * Ls * mb * mt * tA +  4.0 * kmr * Ls * mA * square(mt) * cos(2.0 * (get_bla() - 
		 get_mla())) + 4.0 * kmr * Ls * mA * mb * mt * cos(2.0 * (get_bla() -  get_mla())) - 4.0 * kt * Ls * square(mt) * tA * cos(2.0 * (get_bla() 
		-  get_mla())) - 4.0 * kt * Ls * mb * mt * tA * cos(2.0 * (get_bla() -  get_mla())) - 4.0 * kmr * Ls * square(mt) * M_PI * cos(2.0 
		* (get_bla() -  get_mla())) + 4.0 * kt * Ls * square(mt) * M_PI * cos(2.0 * (get_bla() -  get_mla())) - 4.0 * kmr * Ls * mb * mt * 
		M_PI * cos(2.0 * (get_bla() -  get_mla())) + 4.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * (get_bla() -  get_mla())) - 4.0 * kmr * Lt 
		* mA * square(mt) * cos(get_bra() -  get_mra()) - 4.0 * kmr * Lt * mA * mb * mt * cos(get_bra() -  get_mra()) + 4.0 * kmr * Lt * square(mt) 
		* M_PI * cos(get_bra() -  get_mra()) + 4.0 * kmr * Lt * mb * mt * M_PI * cos(get_bra() -  get_mra()) + 2.0 * bA * kbl * Lt * mb * mt 
		* cos(get_bla() +  2.0 * get_bra() -  get_mra()) - 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mra()) + 
		2.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() +  2.0 * get_bra() -  get_mra()) + 2.0 * bA * kbl * Lt * mb * mt * cos(get_bla() - 
		 2.0 * get_mla() -  get_mra()) - 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) + 2.0 * kml * Lt * 
		mb * mt * M_PI * cos(get_bla() -  2.0 * get_mla() -  get_mra()) + 2.0 * kmr * Lt * mA * square(mt) * cos(2.0 * get_bla() +  get_bra() 
		-  2.0 * get_mla() -  get_mra()) + 2.0 * kmr * Lt * mA * mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) 
		- 2.0 * kmr * Lt * square(mt) * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) - 2.0 * kmr * Lt * mb * mt 
		* M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) - 2.0 * kml * Ls * mA * mb * mt * cos(2.0 * get_bla() + 
		 2.0 * get_bra() -  get_mla() -  get_mra()) + 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() - 
		 get_mra()) + 2.0 * kml * Ls * mb * mt * M_PI * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) - 2.0 * kt * Ls * 
		mb * mt * M_PI * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) + 2.0 * kml * Ls * mA * mb * mt * cos(2.0 * get_bra() 
		+  get_mla() -  get_mra()) - 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bra() +  get_mla() -  get_mra()) - 2.0 * kml * Ls * mb * 
		mt * M_PI * cos(2.0 * get_bra() +  get_mla() -  get_mra()) + 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bra() +  get_mla() -  get_mra()) 
		- 2.0 * bA * kbl * Lt * mb * mt * cos(get_bla() +  get_mra()) + 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() +  get_mra()) - 2.0 * 
		kml * Lt * mb * mt * M_PI * cos(get_bla() +  get_mra()) - 2.0 * bA * kbl * Lt * mb * mt * cos(get_bla() -  2.0 * get_bra() -  2.0 * 
		get_mla() +  get_mra()) + 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) - 2.0 * 
		kml * Lt * mb * mt * M_PI * cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) + 2.0 * kmr * Lt * mA * square(mt) * 
		cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) + 2.0 * kmr * Lt * mA * mb * mt * cos(2.0 * get_bla() -  get_bra() 
		-  2.0 * get_mla() +  get_mra()) - 2.0 * kmr * Lt * square(mt) * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) 
		- 2.0 * kmr * Lt * mb * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) + 2.0 * kml * Ls * mA * mb * 
		mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) - 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() -  get_mla() +  get_mra()) 
		- 2.0 * kml * Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_mla() +  get_mra()) + 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bla() 
		-  get_mla() +  get_mra()) - 2.0 * kml * Ls * mA * mb * mt * cos(get_mla() +  get_mra()) + 2.0 * kt * Ls * mb * mt * tA * cos(get_mla() 
		+  get_mra()) + 2.0 * kml * Ls * mb * mt * M_PI * cos(get_mla() +  get_mra()) - 2.0 * kt * Ls * mb * mt * M_PI * cos(get_mla() +  get_mra()) 
		+ 2.0 * get_bla() * (2.0 * kbl * Lt * cos(get_bla() +  get_mra()) * square(mm) + 2.0 * kml * Lt * cos(get_bla() +  get_mra()) * square(mm) 
		+ 2.0 * kml * Ls * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * square(mm) - 6.0 * kml * Ls * cos(get_mla() +  get_mra()) * square(mm) 
		+ 2.0 * kml * Ls * mb * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * mm + 2.0 * kbl * Lt * mb * cos(get_bla() +  get_mra()) * mm 
		+ 2.0 * kml * Lt * mb * cos(get_bla() +  get_mra()) * mm + 2.0 * kbl * Lt * mt * cos(get_bla() +  get_mra()) * mm + 2.0 * kml * Lt 
		* mt * cos(get_bla() +  get_mra()) * mm + 2.0 * kml * Ls * mb * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm + 2.0 * kml * Ls 
		* mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * mm - 4.0 * kml * Ls * mb * cos(get_mla() +  get_mra()) * mm - 2.0 * kml * Ls 
		* mt * cos(get_mla() +  get_mra()) * mm -  (kbl +  kml) * Lt * mb * mt * cos(get_bla() +  2.0 * get_bra() -  get_mra()) -  (kbl + kml) 
		* Lt * (2.0 * mm * (mm +  mt) +  mb * (2.0 * mm +  mt)) * cos(get_bla() -  2.0 * get_mla() -  get_mra()) - kml * Ls * mb * mt * cos(2.0 
		* get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) + kml * Ls * mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) + kbl 
		* Lt * mb * mt * cos(get_bla() +  get_mra()) + kml * Lt * mb * mt * cos(get_bla() +  get_mra()) + kbl * Lt * mb * mt * cos(get_bla() 
		-  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) + kml * Lt * mb * mt * cos(get_bla() -  2.0 * get_bra() -  2.0 * get_mla() +  get_mra()) 
		+ kml * Ls * mb * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) - kml * Ls * mb * mt * cos(get_mla() +  get_mra())) + 2.0 * kmr 
		* get_bra() * (- 6.0 * Ls * square(mm) -  2.0 * Lt * cos(get_bra() -  get_mra()) * square(mm) + 2.0 * Lt * cos(get_bra() +  2.0 * get_mla() 
		+  get_mra()) * square(mm) -  4.0 * Ls * mb * mm - 8.0 * Ls * mt * mm +  2.0 * Ls * mb * cos(2.0 * (get_bra() +  get_mla())) * mm - 
		2.0 * Lt * mb * cos(get_bra() -  get_mra()) * mm - 6.0 * Lt * mt * cos(get_bra() -  get_mra()) * mm + Lt * mt * cos(2.0 * get_bla() 
		+  get_bra() -  2.0 * get_mla() -  get_mra()) * mm + Lt * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * 
		mm + 2.0 * Lt * mb * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * mm -  2.0 * Ls * square(mt) - 2.0 * Ls * mb * mt + 2.0 * Ls * 
		(mm +  mt) * (mb +  mm +  mt) * cos(2.0 * (get_bla() -  get_mla())) - 2.0 * Lt * square(mt) * cos(get_bra() -  get_mra()) - 2.0 * Lt 
		* mb * mt * cos(get_bra() -  get_mra()) + Lt * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) + Lt 
		* mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) + Lt * square(mt) * cos(2.0 * get_bla() -  get_bra() 
		-  2.0 * get_mla() +  get_mra()) + Lt * mb * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra())) + 4.0 * kt * Ls 
		* square(mt) * get_mla() +  4.0 * kt * Ls * mb * mt * get_mla() - 4.0 * kt * Ls * square(mt) * cos(2.0 * (get_bla() -  get_mla())) 
		* get_mla() - 4.0 * kt * Ls * mb * mt * cos(2.0 * (get_bla() -  get_mla())) * get_mla() + 2.0 * kml * Lt * mb * mt * cos(get_bla() 
		+  2.0 * get_bra() -  get_mra()) * get_mla() + 2.0 * kml * Lt * mb * mt * cos(get_bla() -  2.0 * get_mla() -  get_mra()) * get_mla() 
		+ 2.0 * kml * Ls * mb * mt * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) * get_mla() +  2.0 * kt * Ls * mb * 
		mt * cos(2.0 * get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) * get_mla() - 2.0 * kml * Ls * mb * mt * cos(2.0 * get_bra() 
		+  get_mla() -  get_mra()) * get_mla() - 2.0 * kt * Ls * mb * mt * cos(2.0 * get_bra() +  get_mla() -  get_mra()) * get_mla() - 2.0 
		* kml * Lt * mb * mt * cos(get_bla() +  get_mra()) * get_mla() - 2.0 * kml * Lt * mb * mt * cos(get_bla() -  2.0 * get_bra() -  2.0 
		* get_mla() +  get_mra()) * get_mla() -  2.0 * kml * Ls * mb * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() - 2.0 
		* kt * Ls * mb * mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mla() + 2.0 * kml * Ls * mb * mt * cos(get_mla() +  get_mra()) 
		* get_mla() + 2.0 * kt * Ls * mb * mt * cos(get_mla() +  get_mra()) * get_mla() + 4.0 * kmr * Ls * square(mt) * get_mra() +  4.0 * 
		kt * Ls * square(mt) * get_mra() + 4.0 * kmr * Ls * mb * mt * get_mra() +  4.0 * kt * Ls * mb * mt * get_mra() - 4.0 * kmr * Ls * square(mt) 
		* cos(2.0 * (get_bla() -  get_mla())) * get_mra() - 4.0 * kt * Ls * square(mt) * cos(2.0 * (get_bla() -  get_mla())) * get_mra() - 
		4.0 * kmr * Ls * mb * mt * cos(2.0 * (get_bla() -  get_mla())) * get_mra() - 4.0 * kt * Ls * mb * mt * cos(2.0 * (get_bla() -  get_mla())) 
		* get_mra() + 4.0 * kmr * Lt * square(mt) * cos(get_bra() -  get_mra()) * get_mra() + 4.0 * kmr * Lt * mb * mt * cos(get_bra() -  get_mra()) 
		* get_mra() - 2.0 * kmr * Lt * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() -  2.0 * 
		kmr * Lt * mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() + 2.0 * kt * Ls * mb * mt * cos(2.0 
		* get_bla() +  2.0 * get_bra() -  get_mla() -  get_mra()) * get_mra() -  2.0 * kt * Ls * mb * mt * cos(2.0 * get_bra() +  get_mla() 
		-  get_mra()) * get_mra() - 2.0 * kmr * Lt * square(mt) * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() 
		-  2.0 * kmr * Lt * mb * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() - 2.0 * kt * Ls * mb * 
		mt * cos(2.0 * get_bla() -  get_mla() +  get_mra()) * get_mra() + 2.0 * kt * Ls * mb * mt * cos(get_mla() +  get_mra()) * get_mra() 
		+ 4.0 * kmr * Ls * square(mt) * M_PI -  4.0 * kt * Ls * square(mt) * M_PI + 4.0 * kmr * Ls * mb * mt * M_PI - 4.0 * kt * Ls * mb * 
		mt * M_PI)/(Ls * square(Lt) * (- 4.0 * cos(2.0 * (get_mla() +  get_mra())) * cube(mm) +  4.0 * cube(mm) +  4.0 * mb * square(mm) + 12.0 
		* mt * square(mm) -  4.0 * mb * cos(2.0 * (get_mla() +  get_mra())) * square(mm) + 4.0 * square(mt) * mm +  8.0 * mb * mt * mm +  2.0 
		* mb * square(mt) - 2.0 * mt * (2.0 * mm * (mm +  mt) +  mb * (2.0 * mm +  mt)) * cos(2.0 * (get_bla() -  get_mla())) - 2.0 * mb * 
		mt * (2.0 * mm +  mt) * cos(2.0 * (get_bra() -  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() +  get_bra() -  get_mla() -  get_mra())) 
		+ mb * square(mt) * cos(2.0 * (get_bla() -  get_bra() -  get_mla() +  get_mra()))));
	
	else if (state == RIGHTBOUND) return 0;
	else if (state == BOTHBOUND) return 0;
}

double Dynein::get_dd_bra() {
	if (state == LEFTBOUND) return 
		- (4.0 * kmr * Lt * mA * cube(mm) - 4.0 * kmr * Lt * mA * cos(2.0 * (get_mla() +  get_mra())) * cube(mm) + 4.0 * kmr * Lt * M_PI * cos(2.0 
		* (get_mla() +  get_mra())) * cube(mm) - 4.0 * kmr * Lt * get_mra() * cube(mm) + 4.0 * kmr * Lt * cos(2.0 * (get_mla() +  get_mra())) * 
		get_mra() * cube(mm) - 4.0 * kmr * Lt * M_PI * cube(mm) + 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() +  get_bra()) * square(get_d_bla()) 
		* square(mm) - 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() -  get_bra() -  2.0 * get_mla()) * square(get_d_bla()) * square(mm) 
		- 2.0 * square(Ls) * Lt * mb * mt * sin(get_bla() -  get_bra() +  2.0 * get_mra()) * square(get_d_bla()) * square(mm) + 2.0 * square(Ls) 
		* Lt * mb * mt * sin(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * square(get_d_bla()) * square(mm) + 4.0 * Ls * square(Lt) 
		* mb * mt * sin(get_bra() +  get_mla()) * square(get_d_mla()) * square(mm) + 4.0 * Ls * square(Lt) * mb * mt * sin(get_bra() -  get_mla() 
		-  2.0 * get_mra()) * square(get_d_mla()) * square(mm) + 8.0 * Ls * square(Lt) * mb * mt * sin(get_bra() -  get_mra()) * square(get_d_mra()) 
		* square(mm) +  8.0 * kmr * Lt * mA * mb * square(mm) +  12.0 * kmr * Lt * mA * mt * square(mm) - 4.0 * kmr * Lt * mA * mt * cos(2.0 
		* (get_bla() -  get_mla())) * square(mm) + 4.0 * kmr * Lt * mt * M_PI * cos(2.0 * (get_bla() -  get_mla())) * square(mm) - 4.0 * kml 
		* Ls * mA * mb * cos(get_bra() +  get_mla()) * square(mm) + 4.0 * kt * Ls * mb * tA * cos(get_bra() +  get_mla()) * square(mm) + 4.0 
		* kml * Ls * mb * M_PI * cos(get_bra() +  get_mla()) * square(mm) - 4.0 * kt * Ls * mb * M_PI * cos(get_bra() +  get_mla()) * square(mm) 
		+ 4.0 * kml * Ls * mA * mb * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * square(mm) - 4.0 * kt * Ls * mb * tA * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * square(mm) - 4.0 * kml * Ls * mb * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * square(mm) 
		+ 4.0 * kt * Ls * mb * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * square(mm) + 4.0 * kmr * Ls * mA * mb * cos(get_bra() 
		-  get_mra()) * square(mm) - 4.0 * kt * Ls * mb * tA * cos(get_bra() -  get_mra()) * square(mm) - 4.0 * kmr * Ls * mb * M_PI * cos(get_bra() 
		-  get_mra()) * square(mm) + 4.0 * kt * Ls * mb * M_PI * cos(get_bra() -  get_mra()) * square(mm) - 8.0 * kmr * Lt * mA * mb * cos(2.0 
		* (get_mla() +  get_mra())) * square(mm) + 8.0 * kmr * Lt * mb * M_PI * cos(2.0 * (get_mla() +  get_mra())) * square(mm) - 4.0 * kmr 
		* Ls * mA * mb * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * square(mm) + 4.0 * kt * Ls * mb * tA * cos(get_bra() +  2.0 * get_mla() 
		+  get_mra()) * square(mm) + 4.0 * kmr * Ls * mb * M_PI * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * square(mm) - 4.0 * kt * 
		Ls * mb * M_PI * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * square(mm) + 4.0 * kml * Ls * mb * cos(get_bra() +  get_mla()) * 
		get_mla() * square(mm) + 4.0 * kt * Ls * mb * cos(get_bra() +  get_mla()) * get_mla() * square(mm) - 4.0 * kml * Ls * mb * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * get_mla() * square(mm) - 4.0 * kt * Ls * mb * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mla() 
		* square(mm) - 4.0 * kt * Ls * mb * cos(get_bra() -  get_mra()) * get_mla() * square(mm) + 4.0 * kt * Ls * mb * cos(get_bra() +  2.0 
		* get_mla() +  get_mra()) * get_mla() * square(mm) - 8.0 * kmr * Lt * mb * get_mra() * square(mm) -  12.0 * kmr * Lt * mt * get_mra() 
		* square(mm) + 4.0 * kmr * Lt * mt * cos(2.0 * (get_bla() -  get_mla())) * get_mra() * square(mm) + 4.0 * kt * Ls * mb * cos(get_bra() 
		+  get_mla()) * get_mra() * square(mm) - 4.0 * kt * Ls * mb * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mra() * square(mm) 
		- 4.0 * kmr * Ls * mb * cos(get_bra() -  get_mra()) * get_mra() * square(mm) - 4.0 * kt * Ls * mb * cos(get_bra() -  get_mra()) * get_mra() 
		* square(mm) + 8.0 * kmr * Lt * mb * cos(2.0 * (get_mla() +  get_mra())) * get_mra() * square(mm) + 4.0 * kmr * Ls * mb * cos(get_bra() 
		+  2.0 * get_mla() +  get_mra()) * get_mra() * square(mm) + 4.0 * kt * Ls * mb * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * get_mra() 
		* square(mm) - 8.0 * kmr * Lt * mb * M_PI * square(mm) -  12.0 * kmr * Lt * mt * M_PI * square(mm) + 4.0 * kmr * Lt * mA * square(mb) 
		* mm +  4.0 * kmr * Lt * mA * square(mt) * mm + 2.0 * square(Ls) * Lt * square(mb) * mt * sin(get_bla() +  get_bra()) * square(get_d_bla()) 
		* mm - 2.0 * square(Ls) * Lt * square(mb) * mt * sin(get_bla() -  get_bra() -  2.0 * get_mla()) * square(get_d_bla()) * mm - 2.0 * 
		square(Ls) * Lt * square(mb) * mt * sin(get_bla() -  get_bra() +  2.0 * get_mra()) * square(get_d_bla()) * mm + 2.0 * square(Ls) * 
		Lt * square(mb) * mt * sin(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * square(get_d_bla()) * mm + 4.0 * square(Ls) 
		* Lt * square(mb) * mt * sin(2.0 * (get_bra() -  get_mra())) * square(get_d_bra()) * mm + 4.0 * Ls * square(Lt) * square(mb) * mt * 
		sin(get_bra() +  get_mla()) * square(get_d_mla()) * mm + 4.0 * Ls * square(Lt) * square(mb) * mt * sin(get_bra() -  get_mla() -  2.0 
		* get_mra()) * square(get_d_mla()) * mm + 4.0 * Ls * square(Lt) * mb * square(mt) * sin(get_bra() -  get_mra()) * square(get_d_mra()) 
		* mm + 8.0 * Ls * square(Lt) * square(mb) * mt * sin(get_bra() -  get_mra()) * square(get_d_mra()) * mm - 2.0 * Ls * square(Lt) * mb 
		* square(mt) * sin(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * square(get_d_mra()) * mm + 2.0 * Ls * square(Lt) 
		* mb * square(mt) * sin(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * square(get_d_mra()) * mm +  16.0 * kmr * Lt 
		* mA * mb * mt * mm - 2.0 * bA * kbl * Lt * mb * mt * cos(get_bla() +  get_bra()) * mm + 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() 
		+  get_bra()) * mm - 2.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() +  get_bra()) * mm + 2.0 * bA * kbl * Lt * mb * mt * cos(get_bla() 
		-  get_bra() -  2.0 * get_mla()) * mm - 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) * mm + 2.0 * 
		kml * Lt * mb * mt * M_PI * cos(get_bla() -  get_bra() -  2.0 * get_mla()) * mm - 4.0 * kmr * Lt * mA * square(mt) * cos(2.0 * (get_bla() 
		-  get_mla())) * mm - 8.0 * kmr * Lt * mA * mb * mt * cos(2.0 * (get_bla() -  get_mla())) * mm + 4.0 * kmr * Lt * square(mt) * M_PI 
		* cos(2.0 * (get_bla() -  get_mla())) * mm + 8.0 * kmr * Lt * mb * mt * M_PI * cos(2.0 * (get_bla() -  get_mla())) * mm + 2.0 * kml 
		* Ls * mA * mb * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * mm - 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() +  get_bra() 
		-  get_mla()) * mm - 2.0 * kml * Ls * mb * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * mm + 2.0 * kt * Ls * mb * mt 
		* M_PI * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * mm - 4.0 * kml * Ls * mA * square(mb) * cos(get_bra() +  get_mla()) * mm 
		- 2.0 * kml * Ls * mA * mb * mt * cos(get_bra() +  get_mla()) * mm + 4.0 * kt * Ls * square(mb) * tA * cos(get_bra() +  get_mla()) 
		* mm + 2.0 * kt * Ls * mb * mt * tA * cos(get_bra() +  get_mla()) * mm + 4.0 * kml * Ls * square(mb) * M_PI * cos(get_bra() +  get_mla()) 
		* mm - 4.0 * kt * Ls * square(mb) * M_PI * cos(get_bra() +  get_mla()) * mm + 2.0 * kml * Ls * mb * mt * M_PI * cos(get_bra() +  get_mla()) 
		* mm - 2.0 * kt * Ls * mb * mt * M_PI * cos(get_bra() +  get_mla()) * mm + 4.0 * kml * Ls * mA * square(mb) * cos(get_bra() -  get_mla() 
		-  2.0 * get_mra()) * mm + 2.0 * kml * Ls * mA * mb * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm - 4.0 * kt * Ls * square(mb) 
		* tA * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm - 2.0 * kt * Ls * mb * mt * tA * cos(get_bra() -  get_mla() -  2.0 * get_mra()) 
		* mm - 4.0 * kml * Ls * square(mb) * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm + 4.0 * kt * Ls * square(mb) * M_PI 
		* cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm - 2.0 * kml * Ls * mb * mt * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) 
		* mm + 2.0 * kt * Ls * mb * mt * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * mm + 4.0 * kmr * Ls * mA * square(mb) * cos(get_bra() 
		-  get_mra()) * mm + 12.0 * kmr * Ls * mA * mb * mt * cos(get_bra() -  get_mra()) * mm - 4.0 * kt * Ls * square(mb) * tA * cos(get_bra() 
		-  get_mra()) * mm - 12.0 * kt * Ls * mb * mt * tA * cos(get_bra() -  get_mra()) * mm - 4.0 * kmr * Ls * square(mb) * M_PI * cos(get_bra() 
		-  get_mra()) * mm + 4.0 * kt * Ls * square(mb) * M_PI * cos(get_bra() -  get_mra()) * mm - 12.0 * kmr * Ls * mb * mt * M_PI * cos(get_bra() 
		-  get_mra()) * mm + 12.0 * kt * Ls * mb * mt * M_PI * cos(get_bra() -  get_mra()) * mm - 2.0 * kmr * Ls * mA * mb * mt * cos(2.0 * 
		get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * mm + 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() +  get_bra() -  2.0 
		* get_mla() -  get_mra()) * mm + 2.0 * kmr * Ls * mb * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) 
		* mm - 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * mm - 2.0 * kmr * Ls * mA 
		* mb * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * mm + 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() 
		-  get_bra() -  2.0 * get_mla() +  get_mra()) * mm + 2.0 * kmr * Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() 
		+  get_mra()) * mm - 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * mm - 4.0 
		* kmr * Lt * mA * square(mb) * cos(2.0 * (get_mla() +  get_mra())) * mm + 4.0 * kmr * Lt * square(mb) * M_PI * cos(2.0 * (get_mla() 
		+  get_mra())) * mm - 4.0 * kmr * Ls * mA * square(mb) * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * mm + 4.0 * kt * Ls * square(mb) 
		* tA * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * mm + 4.0 * kmr * Ls * square(mb) * M_PI * cos(get_bra() +  2.0 * get_mla() 
		+  get_mra()) * mm - 4.0 * kt * Ls * square(mb) * M_PI * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * mm + 2.0 * bA * kbl * Lt 
		* mb * mt * cos(get_bla() -  get_bra() +  2.0 * get_mra()) * mm - 2.0 * kml * Lt * mA * mb * mt * cos(get_bla() -  get_bra() +  2.0 
		* get_mra()) * mm + 2.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() -  get_bra() +  2.0 * get_mra()) * mm - 2.0 * kml * Ls * mA * mb 
		* mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * mm + 2.0 * kt * Ls * mb * mt * tA * cos(2.0 * get_bla() 
		-  get_bra() -  get_mla() +  2.0 * get_mra()) * mm + 2.0 * kml * Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  get_mla() 
		+  2.0 * get_mra()) * mm - 2.0 * kt * Ls * mb * mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * mm 
		- 2.0 * bA * kbl * Lt * mb * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * mm + 2.0 * kml * Lt * mA * mb * mt 
		* cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * mm - 2.0 * kml * Lt * mb * mt * M_PI * cos(get_bla() +  get_bra() 
		-  2.0 * (get_mla() +  get_mra())) * mm - 2.0 * kml * Lt * mb * mt * cos(get_bla() +  get_bra()) * get_mla() * mm + 2.0 * kml * Lt 
		* mb * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) * get_mla() * mm - 2.0 * kml * Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() 
		-  get_mla()) * get_mla() * mm - 2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * get_mla() * mm + 4.0 * 
		kml * Ls * square(mb) * cos(get_bra() +  get_mla()) * get_mla() * mm + 4.0 * kt * Ls * square(mb) * cos(get_bra() +  get_mla()) * get_mla() 
		* mm + 2.0 * kml * Ls * mb * mt * cos(get_bra() +  get_mla()) * get_mla() * mm + 2.0 * kt * Ls * mb * mt * cos(get_bra() +  get_mla()) 
		* get_mla() * mm - 4.0 * kml * Ls * square(mb) * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mla() * mm - 4.0 * kt * Ls * 
		square(mb) * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mla() * mm - 2.0 * kml * Ls * mb * mt * cos(get_bra() -  get_mla() 
		-  2.0 * get_mra()) * get_mla() * mm - 2.0 * kt * Ls * mb * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mla() * mm - 
		4.0 * kt * Ls * square(mb) * cos(get_bra() -  get_mra()) * get_mla() * mm - 12.0 * kt * Ls * mb * mt * cos(get_bra() -  get_mra()) 
		* get_mla() * mm + 2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mla() * mm +  
		2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mla() * mm + 4.0 * kt * Ls * square(mb) 
		* cos(get_bra() +  2.0 * get_mla() +  get_mra()) * get_mla() * mm + 2.0 * kml * Lt * mb * mt * cos(get_bla() -  get_bra() +  2.0 * 
		get_mra()) * get_mla() * mm + 2.0 * kml * Ls * mb * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * get_mla() 
		* mm +  2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * get_mla() * mm - 2.0 * kml * 
		Lt * mb * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * get_mla() * mm -  4.0 * kmr * Lt * square(mb) * get_mra() 
		* mm -  4.0 * kmr * Lt * square(mt) * get_mra() * mm - 16.0 * kmr * Lt * mb * mt * get_mra() * mm + 4.0 * kmr * Lt * square(mt) * cos(2.0 
		* (get_bla() -  get_mla())) * get_mra() * mm + 8.0 * kmr * Lt * mb * mt * cos(2.0 * (get_bla() -  get_mla())) * get_mra() * mm - 2.0 
		* kt * Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * get_mra() * mm + 4.0 * kt * Ls * square(mb) * cos(get_bra() 
		+  get_mla()) * get_mra() * mm + 2.0 * kt * Ls * mb * mt * cos(get_bra() +  get_mla()) * get_mra() * mm - 4.0 * kt * Ls * square(mb) 
		* cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mra() * mm - 2.0 * kt * Ls * mb * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) 
		* get_mra() * mm - 4.0 * kmr * Ls * square(mb) * cos(get_bra() -  get_mra()) * get_mra() * mm - 4.0 * kt * Ls * square(mb) * cos(get_bra() 
		-  get_mra()) * get_mra() * mm - 12.0 * kmr * Ls * mb * mt * cos(get_bra() -  get_mra()) * get_mra() * mm - 12.0 * kt * Ls * mb * mt 
		* cos(get_bra() -  get_mra()) * get_mra() * mm + 2.0 * kmr * Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  
		get_mra()) * get_mra() * mm +  2.0 * kt * Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() 
		* mm + 2.0 * kmr * Ls * mb * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() * mm +  2.0 * kt * 
		Ls * mb * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() * mm + 4.0 * kmr * Lt * square(mb) * cos(2.0 
		* (get_mla() +  get_mra())) * get_mra() * mm + 4.0 * kmr * Ls * square(mb) * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * get_mra() 
		* mm + 4.0 * kt * Ls * square(mb) * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * get_mra() * mm + 2.0 * kt * Ls * mb * mt * cos(2.0 
		* get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * get_mra() * mm -  4.0 * kmr * Lt * square(mb) * M_PI * mm -  4.0 * kmr 
		* Lt * square(mt) * M_PI * mm - 16.0 * kmr * Lt * mb * mt * M_PI * mm +  4.0 * kmr * Lt * mA * mb * square(mt) + 2.0 * square(Ls) * 
		Lt * square(mb) * square(mt) * sin(2.0 * (get_bra() -  get_mra())) * square(get_d_bra()) - square(Ls) * Lt * square(mb) * square(mt) 
		* sin(2.0 * (get_bla() +  get_bra() -  get_mla() -  get_mra())) * square(get_d_bra()) +  square(Ls) * Lt * square(mb) * square(mt) 
		* sin(2.0 * (get_bla() -  get_bra() -  get_mla() +  get_mra())) * square(get_d_bra()) + 4.0 * Ls * square(Lt) * square(mb) * square(mt) 
		* sin(get_bra() -  get_mra()) * square(get_d_mra()) -  2.0 * Ls * square(Lt) * square(mb) * square(mt) * sin(2.0 * get_bla() +  get_bra() 
		-  2.0 * get_mla() -  get_mra()) * square(get_d_mra()) +  2.0 * Ls * square(Lt) * square(mb) * square(mt) * sin(2.0 * get_bla() -  
		get_bra() -  2.0 * get_mla() +  get_mra()) * square(get_d_mra()) +  4.0 * kmr * Lt * mA * square(mb) * mt - 2.0 * bA * kbl * Lt * square(mb) 
		* mt * cos(get_bla() +  get_bra()) + 2.0 * kml * Lt * mA * square(mb) * mt * cos(get_bla() +  get_bra()) - 2.0 * kml * Lt * square(mb) 
		* mt * M_PI * cos(get_bla() +  get_bra()) + 2.0 * bA * kbl * Lt * square(mb) * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) 
		- 2.0 * kml * Lt * mA * square(mb) * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) + 2.0 * kml * Lt * square(mb) * mt * M_PI 
		* cos(get_bla() -  get_bra() -  2.0 * get_mla()) - 4.0 * kmr * Lt * mA * mb * square(mt) * cos(2.0 * (get_bla() -  get_mla())) - 4.0 
		* kmr * Lt * mA * square(mb) * mt * cos(2.0 * (get_bla() -  get_mla())) + 4.0 * kmr * Lt * mb * square(mt) * M_PI * cos(2.0 * (get_bla() 
		-  get_mla())) + 4.0 * kmr * Lt * square(mb) * mt * M_PI * cos(2.0 * (get_bla() -  get_mla())) + 2.0 * kml * Ls * mA * square(mb) * 
		mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) - 2.0 * kt * Ls * square(mb) * mt * tA * cos(2.0 * get_bla() +  get_bra() -  get_mla()) 
		- 2.0 * kml * Ls * square(mb) * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  get_mla()) + 2.0 * kt * Ls * square(mb) * mt * M_PI 
		* cos(2.0 * get_bla() +  get_bra() -  get_mla()) - 2.0 * kml * Ls * mA * square(mb) * mt * cos(get_bra() +  get_mla()) + 2.0 * kt * 
		Ls * square(mb) * mt * tA * cos(get_bra() +  get_mla()) + 2.0 * kml * Ls * square(mb) * mt * M_PI * cos(get_bra() +  get_mla()) - 2.0 
		* kt * Ls * square(mb) * mt * M_PI * cos(get_bra() +  get_mla()) + 2.0 * kml * Ls * mA * square(mb) * mt * cos(get_bra() -  get_mla() 
		-  2.0 * get_mra()) - 2.0 * kt * Ls * square(mb) * mt * tA * cos(get_bra() -  get_mla() -  2.0 * get_mra()) - 2.0 * kml * Ls * square(mb) 
		* mt * M_PI * cos(get_bra() -  get_mla() -  2.0 * get_mra()) + 2.0 * kt * Ls * square(mb) * mt * M_PI * cos(get_bra() -  get_mla() 
		-  2.0 * get_mra()) + 4.0 * kmr * Ls * mA * mb * square(mt) * cos(get_bra() -  get_mra()) + 4.0 * kmr * Ls * mA * square(mb) * mt * 
		cos(get_bra() -  get_mra()) - 4.0 * kt * Ls * mb * square(mt) * tA * cos(get_bra() -  get_mra()) - 4.0 * kt * Ls * square(mb) * mt 
		* tA * cos(get_bra() -  get_mra()) - 4.0 * kmr * Ls * mb * square(mt) * M_PI * cos(get_bra() -  get_mra()) + 4.0 * kt * Ls * mb * square(mt) 
		* M_PI * cos(get_bra() -  get_mra()) - 4.0 * kmr * Ls * square(mb) * mt * M_PI * cos(get_bra() -  get_mra()) + 4.0 * kt * Ls * square(mb) 
		* mt * M_PI * cos(get_bra() -  get_mra()) - 2.0 * kmr * Ls * mA * mb * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() 
		-  get_mra()) - 2.0 * kmr * Ls * mA * square(mb) * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) + 2.0 * kt 
		* Ls * mb * square(mt) * tA * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) + 2.0 * kt * Ls * square(mb) * mt * 
		tA * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) + 2.0 * kmr * Ls * mb * square(mt) * M_PI * cos(2.0 * get_bla() 
		+  get_bra() -  2.0 * get_mla() -  get_mra()) - 2.0 * kt * Ls * mb * square(mt) * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * 
		get_mla() -  get_mra()) + 2.0 * kmr * Ls * square(mb) * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) 
		- 2.0 * kt * Ls * square(mb) * mt * M_PI * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) - 2.0 * kmr * Ls * mA 
		* mb * square(mt) * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) - 2.0 * kmr * Ls * mA * square(mb) * mt * cos(2.0 
		* get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) + 2.0 * kt * Ls * mb * square(mt) * tA * cos(2.0 * get_bla() -  get_bra() 
		-  2.0 * get_mla() +  get_mra()) + 2.0 * kt * Ls * square(mb) * mt * tA * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) 
		+ 2.0 * kmr * Ls * mb * square(mt) * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) - 2.0 * kt * Ls * mb 
		* square(mt) * M_PI * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) + 2.0 * kmr * Ls * square(mb) * mt * M_PI * 
		cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) - 2.0 * kt * Ls * square(mb) * mt * M_PI * cos(2.0 * get_bla() - 
		 get_bra() -  2.0 * get_mla() +  get_mra()) - 2.0 * kmr * get_bra() * (- 2.0 * Lt * cube(mm) +  2.0 * Lt * cos(2.0 * (get_mla() +  get_mra())) 
		* cube(mm) - 4.0 * Lt * mb * square(mm) -  6.0 * Lt * mt * square(mm) + 4.0 * Lt * mb * cos(2.0 * (get_mla() +  get_mra())) * square(mm) 
		+ 2.0 * Ls * mb * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * square(mm) - 2.0 * Lt * square(mb) * mm -  2.0 * Lt * square(mt) 
		* mm -  8.0 * Lt * mb * mt * mm + Ls * mb * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * mm + Ls * mb * 
		mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * mm + 2.0 * Lt * square(mb) * cos(2.0 * (get_mla() +  get_mra())) 
		* mm + 2.0 * Ls * square(mb) * cos(get_bra() +  2.0 * get_mla() +  get_mra()) * mm - 2.0 * Lt * mb * square(mt) -  2.0 * Lt * square(mb) 
		* mt + 2.0 * Lt * (mb +  mm) * mt * (mb +  mm +  mt) * cos(2.0 * (get_bla() -  get_mla())) - 2.0 * Ls * mb * (square(mm) +  3.0 * mt 
		* mm +  square(mt) +  mb * (mm +  mt)) * cos(get_bra() -  get_mra()) + Ls * mb * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 
		* get_mla() -  get_mra()) + Ls * square(mb) * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) + Ls * mb * square(mt) 
		* cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) + Ls * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 
		* get_mla() +  get_mra())) + 2.0 * bA * kbl * Lt * square(mb) * mt * cos(get_bla() -  get_bra() +  2.0 * get_mra()) - 2.0 * kml * Lt 
		* mA * square(mb) * mt * cos(get_bla() -  get_bra() +  2.0 * get_mra()) + 2.0 * kml * Lt * square(mb) * mt * M_PI * cos(get_bla() - 
		 get_bra() +  2.0 * get_mra()) - 2.0 * kml * Ls * mA * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) 
		+ 2.0 * kt * Ls * square(mb) * mt * tA * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) + 2.0 * kml * Ls * square(mb) 
		* mt * M_PI * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) - 2.0 * kt * Ls * square(mb) * mt * M_PI * cos(2.0 
		* get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) - 2.0 * bA * kbl * Lt * square(mb) * mt * cos(get_bla() +  get_bra() -  2.0 
		* (get_mla() +  get_mra())) + 2.0 * kml * Lt * mA * square(mb) * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) 
		- 2.0 * kml * Lt * square(mb) * mt * M_PI * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) - 2.0 * kml * Lt * square(mb) 
		* mt * cos(get_bla() +  get_bra()) * get_mla() + 2.0 * kml * Lt * square(mb) * mt * cos(get_bla() -  get_bra() -  2.0 * get_mla()) 
		* get_mla() - 2.0 * kml * Ls * square(mb) * mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * get_mla() - 2.0 * kt * Ls * square(mb) 
		* mt * cos(2.0 * get_bla() +  get_bra() -  get_mla()) * get_mla() + 2.0 * kml * Ls * square(mb) * mt * cos(get_bra() +  get_mla()) 
		* get_mla() + 2.0 * kt * Ls * square(mb) * mt * cos(get_bra() +  get_mla()) * get_mla() - 2.0 * kml * Ls * square(mb) * mt * cos(get_bra() 
		-  get_mla() -  2.0 * get_mra()) * get_mla() - 2.0 * kt * Ls * square(mb) * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mla() 
		- 4.0 * kt * Ls * mb * square(mt) * cos(get_bra() -  get_mra()) * get_mla() - 4.0 * kt * Ls * square(mb) * mt * cos(get_bra() -  get_mra()) 
		* get_mla() + 2.0 * kt * Ls * mb * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mla() +  2.0 
		* kt * Ls * square(mb) * mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mla() + 2.0 * kt * Ls * mb * 
		square(mt) * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mla() +  2.0 * kt * Ls * square(mb) * mt * cos(2.0 
		* get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mla() + 2.0 * kml * Lt * square(mb) * mt * cos(get_bla() -  get_bra() 
		+  2.0 * get_mra()) * get_mla() + 2.0 * kml * Ls * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) 
		* get_mla() + 2.0 * kt * Ls * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * get_mla() -  2.0 
		* kml * Lt * square(mb) * mt * cos(get_bla() +  get_bra() -  2.0 * (get_mla() +  get_mra())) * get_mla() - 4.0 * kmr * Lt * mb * square(mt) 
		* get_mra() -  4.0 * kmr * Lt * square(mb) * mt * get_mra() + 4.0 * kmr * Lt * mb * square(mt) * cos(2.0 * (get_bla() -  get_mla())) 
		* get_mra() + 4.0 * kmr * Lt * square(mb) * mt * cos(2.0 * (get_bla() -  get_mla())) * get_mra() - 2.0 * kt * Ls * square(mb) * mt 
		* cos(2.0 * get_bla() +  get_bra() -  get_mla()) * get_mra() + 2.0 * kt * Ls * square(mb) * mt * cos(get_bra() +  get_mla()) * get_mra() 
		- 2.0 * kt * Ls * square(mb) * mt * cos(get_bra() -  get_mla() -  2.0 * get_mra()) * get_mra() - 4.0 * kmr * Ls * mb * square(mt) * 
		cos(get_bra() -  get_mra()) * get_mra() - 4.0 * kt * Ls * mb * square(mt) * cos(get_bra() -  get_mra()) * get_mra() - 4.0 * kmr * Ls 
		* square(mb) * mt * cos(get_bra() -  get_mra()) * get_mra() - 4.0 * kt * Ls * square(mb) * mt * cos(get_bra() -  get_mra()) * get_mra() 
		+ 2.0 * kmr * Ls * mb * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() + 2.0 * kt * Ls 
		* mb * square(mt) * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() +  2.0 * kmr * Ls * square(mb) * 
		mt * cos(2.0 * get_bla() +  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() + 2.0 * kt * Ls * square(mb) * mt * cos(2.0 * get_bla() 
		+  get_bra() -  2.0 * get_mla() -  get_mra()) * get_mra() +  2.0 * kmr * Ls * mb * square(mt) * cos(2.0 * get_bla() -  get_bra() - 
		 2.0 * get_mla() +  get_mra()) * get_mra() + 2.0 * kt * Ls * mb * square(mt) * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() 
		+  get_mra()) * get_mra() +  2.0 * kmr * Ls * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * 
		get_mra() + 2.0 * kt * Ls * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  2.0 * get_mla() +  get_mra()) * get_mra() +  2.0 
		* kt * Ls * square(mb) * mt * cos(2.0 * get_bla() -  get_bra() -  get_mla() +  2.0 * get_mra()) * get_mra() - 4.0 * mb * (mb +  mm) 
		* get_bla() * sin(get_bra() - get_mra()) * ((kbl +  kml) * Lt * mt * sin(get_bla() -  2.0 * get_mla() -  get_mra()) +  (kbl +  kml) 
		* Lt * mt * sin(get_bla() +  get_mra()) + kml * Ls * (mt * sin(2.0 * get_bla() -  get_mla() +  get_mra()) -  (2.0 * mm +  mt) * sin(get_mla() 
		+  get_mra()))) -  4.0 * kmr * Lt * mb * square(mt) * M_PI - 4.0 * kmr * Lt * square(mb) * mt * M_PI)/(square(Ls) * Lt * mb * (- 4.0 
		* cos(2.0 * (get_mla() +  get_mra())) * cube(mm) +  4.0 * cube(mm) +  4.0 * mb * square(mm) + 12.0 * mt * square(mm) -  4.0 * mb * cos(2.0 
		* (get_mla() +  get_mra())) * square(mm) + 4.0 * square(mt) * mm +  8.0 * mb * mt * mm +  2.0 * mb * square(mt) - 2.0 * mt * (2.0 * 
		mm * (mm +  mt) +  mb * (2.0 * mm +  mt)) * cos(2.0 * (get_bla() -  get_mla())) - 2.0 * mb * mt * (2.0 * mm +  mt) * cos(2.0 * (get_bra() 
		-  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() +  get_bra() -  get_mla() -  get_mra())) + mb * square(mt) * cos(2.0 * (get_bla() 
		-  get_bra() -  get_mla() +  get_mra()))));
	
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
	return (1.0/2) * kbl * pow(get_bla() - bA, 2) + (1.0/2) * kml * pow(M_PI - get_bla() + get_mla() - mA, 2) + (1.0/2) * kt * 
		pow(M_PI - get_mla() - get_mra() - tA, 2) + (1.0/2) * kmr * pow(M_PI - get_bra() + get_mra() - mA, 2);
}

double Dynein::get_KE() {
	return 1.0/2.0 * mb * (square(- (Ls * get_d_bla() * sin(get_bla())) - Ls * get_d_bra() * sin(get_bra()) - Lt * get_d_mla() * sin(get_mla()) 
		- Lt * get_d_mra() * sin(get_mra())) +  square(Ls * get_d_bla() * cos(get_bla()) -  Ls * get_d_bra() * cos(get_bra()) + Lt * get_d_mla() 
		* cos(get_mla()) - Lt * get_d_mra() * cos(get_mra()))) + 1.0/2.0 * mm * (square(Ls) * square(get_d_bla()) * square(sin(get_bla())) 
		+ square(Ls) * square(get_d_bla()) * square(cos(get_bla()))) + 1.0/2.0 * mm * (square(- (Ls * get_d_bla() * sin(get_bla())) - Lt * 
		get_d_mla() * sin(get_mla()) - Lt * get_d_mra() * sin(get_mra())) +  square(Ls * get_d_bla() * cos(get_bla()) +  Lt * get_d_mla() * 
		cos(get_mla()) - Lt * get_d_mra() * cos(get_mra()))) + 1.0/2.0 * mt * (square(- (Ls * get_d_bla() * sin(get_bla())) - Lt * get_d_mla() 
		* sin(get_mla())) +  square(Ls * get_d_bla() * cos(get_bla()) +  Lt * get_d_mla() * cos(get_mla())));
	
	//Rewrite with new d_blx, etc.
}

void Dynein::log(double t) {
	fprintf(data_file, "%.6f\t%12.6f\t%12.6f\t%.3f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%+.5f\t%d\n", 
	get_KE(), get_PE(), get_KE() + get_PE(), t, get_blx(), get_bly(), get_mlx(), get_mly(), get_tx(), get_ty(), get_mrx(), 
	get_mry(), get_brx(), get_bry(), get_state());
}

Dynein* initProtein(Dynein* dyn) {
	dyn->set_blx(0);
	dyn->set_bly(0);
	
	dyn->set_d_bla(0);
	dyn->set_d_mla(0);
	dyn->set_d_mra(0);
	dyn->set_d_bra(0);
	
	dyn->set_bla(bla_init);
	dyn->set_mla(mla_init);
	dyn->set_mra(mra_init);
	dyn->set_bra(bra_init);
}

Dynein* simulateProtein(Dynein* dyn, double dt, double tf) {
	double t = 0;
	//Change floats to const double, put temp declarations inside while loop
	double temp_bla;
	double temp_mla;
	double temp_mra;
	double temp_bra;
	
	double temp_d_bla = 0;
	double temp_d_mla = 0;
	double temp_d_mra = 0;
	double temp_d_bra = 0;
	
	double temp_dd_bla = 0;
	double temp_dd_mla = 0;
	double temp_dd_mra = 0;
	double temp_dd_bra = 0;
	
	while( t < tf ) {
		temp_d_bla = dyn->get_d_bla() + dyn->get_dd_bla() * dt;
		temp_d_mla = dyn->get_d_mla() + dyn->get_dd_mla() * dt;
		temp_d_mra = dyn->get_d_mra() + dyn->get_dd_mra() * dt;
		temp_d_bra = dyn->get_d_bra() + dyn->get_dd_bra() * dt;
		
		dyn->set_d_bla(temp_d_bla);
		dyn->set_d_mla(temp_d_mla);
		dyn->set_d_mra(temp_d_mra);
		dyn->set_d_bra(temp_d_bra);
		
		temp_bla = dyn->get_bla() + dyn->get_d_bla() * dt;
		temp_mla = dyn->get_mla() + dyn->get_d_mla() * dt;
		temp_mra = dyn->get_mra() + dyn->get_d_mra() * dt;
		temp_bra = dyn->get_bra() + dyn->get_d_bra() * dt;
		
		dyn->set_bla(temp_bla);
		dyn->set_mla(temp_mla);
		dyn->set_mra(temp_mra);
		dyn->set_bra(temp_bra);
		
		if (dyn->get_state() == LEFTBOUND) {
			if (dyn->get_brx() <= -0.1)
				dyn->set_state(BOTHBOUND);
		}
		
		if (dyn->get_state() == RIGHTBOUND) {
		}
		
		if (dyn->get_state() == BOTHBOUND) {
			
		}
																																																																																																			
		dyn->log(t);
		
		t += dt;
	}
	
}

/* *********************************** UTILITY ****************************************** */

double dist(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1) * range;
}

void resetLog(Dynein* dyn) {
	data_file = fopen("data.txt", "w");
	config_file = fopen("config.txt", "w");
	
	fprintf(config_file, "%+.3f\t%+.3f\t%d\n", inctime, runtime, (int) dyn->get_state());
	fprintf(data_file,
		"#KE\t\t\t\tPE\t\t\t\tEnergy\t\tt\t\tblX\t\t\tblY\t\t\tmlX\t\t\tmlY\t\t\ttX\t\t\ttY\t\t\tmrX\t\t\tmrY\t\t\tbrX\t\t\tbrY\t\t\tS\n");
}

double square(double num) {
	return num * num;
}

double cube(double num) {
	return num * num * num;
}

/* *********************************** MAIN ****************************************** */

int main() {
	Dynein* dyn = (Dynein*) malloc(sizeof(Dynein));
	resetLog(dyn);
	initProtein(dyn);
	simulateProtein(dyn, inctime, runtime);
	free(dyn);
	dyn = NULL;
	return 0;
}
