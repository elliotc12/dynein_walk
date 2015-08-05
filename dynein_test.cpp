#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"


int main() {
	
	float bla_t = (108.0 / 180) * M_PI;
	float mla_t = (36.0 / 180) * M_PI;
	float mra_t = (36.0 / 180) * M_PI;
	float bra_t = (108.0 / 180) * M_PI;
	
	Dynein* dyn = new Dynein(bla_t, mla_t, mra_t, bra_t);
	
	printf("one: %f\n", (dyn->get_d_mrx() + lt*sin(dyn->get_bra() - M_PI)*dyn->get_d_bra() - (1/g)*dyn->get_f_brx() - dyn->get_r_brx() ) * 1/(dyn->get_brx() - dyn->get_mrx()));
	printf("two: %f\n", (dyn->get_d_mry() - lt*cos(dyn->get_bra() - M_PI)*dyn->get_d_bra() - (1/g)*dyn->get_f_bry() - dyn->get_r_bry() ) * 1/(dyn->get_bry() - dyn->get_mry()));
	
	exit(EXIT_SUCCESS);
}













