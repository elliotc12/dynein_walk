#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"

	/*
	 * For every timestep, call update_protein to update internal velocities.
	 * Then do Euler's method to update internal coordinates and log output.
	 * update_protein must be called after every set_x command to update
	 * internal velocities.
	 */

void simulateProtein(Dynein* dyn, double dt, double tf) {
	double t = 0;
	
	double temp_bla;
	double temp_mla;
	double temp_mra;
	double temp_bra;
	
	while( t < tf ) {
		
		dyn->update_protein();
		
		temp_bla = dyn->get_bla() + dyn->get_d_bla() * dt;
		temp_mla = dyn->get_mla() + dyn->get_d_mla() * dt;
		temp_mra = dyn->get_mra() + dyn->get_d_mra() * dt;
		temp_bra = dyn->get_bra() + dyn->get_d_bra() * dt;
		
		dyn->set_bla(temp_bla);
		dyn->set_mla(temp_mla);
		dyn->set_mra(temp_mra);
		dyn->set_bra(temp_bra);
		
		dyn->log(t);
		
		t += dt;
	}
}


/* *********************************** MAIN ****************************************** */

int main(int argvc, char **argv) {
	Dynein* dyn = new Dynein(bla_init, mla_init, mra_init, bra_init);
	resetLog(dyn);
	simulateProtein(dyn, inctime, runtime);
	free(dyn);
	dyn = NULL;
	return 0;
}
