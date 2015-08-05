#include <stdio.h>
#include <stdlib.h>

#include "dynein_struct.h"

	/*
	 * For every timestep, call update_velocities to update internal velocities.
	 * Then do Euler's method to update internal coordinates and log output.
	 * update_velocities must be called after every set_x command to update
	 * internal velocities.
	 */

void simulateProtein(Dynein* dyn, double dt, double tf) {
	double t = 0;
	
	double temp_bla;
	double temp_mla;
	double temp_mra;
	double temp_bra;
	
	while( t < tf ) {
		
		dyn->update_velocities();
		
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
	
	forces f;
	f.r_blx = r_blx_init;     f.r_bly = r_bly_init;
	f.r_mlx = r_mlx_init;     f.r_mly = r_mly_init;
	f.r_tx  =  r_tx_init;	  f.r_ty  =  r_ty_init;
	f.r_mrx = r_mrx_init;     f.r_mry = r_mry_init;
	f.r_brx = r_brx_init;     f.r_bry = r_bry_init;
	
	f.f_blx = f_blx_init;     f.f_bly = f_bly_init;
	f.f_mlx = f_mlx_init;     f.f_mly = f_mly_init;
	f.f_tx  =  f_tx_init;     f.f_ty  =  f_ty_init;
	f.f_mrx = f_mrx_init;     f.f_mry = f_mry_init;
	f.f_brx = f_brx_init;     f.f_bry = f_bry_init;
	
	Dynein* dyn = new Dynein(bla_init, mla_init, mra_init, bra_init, f);
	
	resetLog(dyn);
	simulateProtein(dyn, inctime, runtime);
	free(dyn);
	dyn = NULL;
	return 0;
}
