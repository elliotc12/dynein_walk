#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

//Display output of generated dynein_motion_functions.cpp file for various conditions
//For use in testing accuracy of replace.py translation 

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

int main() {
	const double lt = 10.0;
	const double ls = 10.0;

	const double kt  = 1.0;
	const double kml = 1.0;
	const double kmr = 1.0;
	const double kbl = 1.0;
	const double kbr = 1.0;
		
	const double mb = 1.0;
	const double mm = 1.0;
	const double mt = 1.0;

	const double ba = (108.0 / 180) * M_PI;
	const double ma = (108.0 / 180) * M_PI;
	const double ta = (108.0 / 180) * M_PI;

	const double bla_value = (108.0 / 180) * M_PI;
	const double mla_value = (36.0 / 180) * M_PI;
	const double mra_value = (36.0 / 180) * M_PI;
	const double bra_value = (100.0 / 180) * M_PI;

	const double bla_d_value = 0;
	const double mla_d_value = 0;
	const double mra_d_value = 0;
	const double bra_d_value = 0;
	
	Dynein* dyn = (Dynein*) malloc(sizeof(Dynein));
	initProtein(dyn);
	
	printf("Ke: %f\n", 	   dyn->get_KE());
	printf("dd_bla: %E\n", dyn->get_dd_bla());
	printf("dd_mla: %E\n", dyn->get_dd_mla());
	printf("dd_mra: %E\n", dyn->get_dd_mra());
	printf("dd_bra: %E\n", dyn->get_dd_bra());
	
	free(dyn);
	dyn = NULL;
	
}




//Add functionality into solve.py to output the outputs of its dd functions to this program in a varg which will compare them with the replace.py'd function outputs for error checking.
