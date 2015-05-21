#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

//Display output of generated dynein_motion_functions.cpp file for various conditions
//For use in testing accuracy of replace.py translation

     double fblx = 0;	 double fbly = 0;
	 double fmlx = 0;	 double fmly = 0;
	 double ftx = 0;	 double fty = 0;
	 double fmrx = 0;	 double fmry = 0;
	 double fbrx = 0;	 double fbry = 0;

Dynein* initProtein(Dynein* dyn, double bla_value, double mla_value, double mra_value, double bra_value) {
	dyn->set_blx(0);
	dyn->set_bly(0);
	
	dyn->set_d_bla(0);
	dyn->set_d_mla(0);
	dyn->set_d_mra(0);
	dyn->set_d_bra(0);
	
	dyn->set_bla(bla_value);
	dyn->set_mla(mla_value);
	dyn->set_mra(mra_value);
	dyn->set_bra(bra_value);
}

int main() {
	
	const double lt = 10.0;
	const double ls = 10.0;

	const double bla_value = (108.0 / 180) * M_PI;
	const double mla_value = (36.0 / 180) * M_PI;
	const double mra_value = (36.0 / 180) * M_PI;
	const double bra_value = (108.0 / 180) * M_PI;
	
	Dynein* dyn = (Dynein*) malloc(sizeof(Dynein));
	initProtein(dyn, bla_value, mla_value, mra_value, bra_value);
	
	dyn->set_state(LEFTBOUND);
	printf("\n\nC Leftbound Accelerations:\n\n");
	printf("dd_bla: %E\n", dyn->get_dd_bla());
	printf("dd_mla: %E\n", dyn->get_dd_mla());
	printf("dd_mra: %E\n", dyn->get_dd_mra());
	printf("dd_bra: %E\n", dyn->get_dd_bra());
	
	free(dyn);
	dyn = NULL;
	
}


