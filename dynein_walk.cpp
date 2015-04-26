#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"

double fblx = 0;		double fbly = 0;
double fmlx = 0;		double fmly = 0;
double ftx = 0;			double fty = 0;
double fmrx = 0;		double fmry = 0;
double fbrx = 0;		double fbry = 0;

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
