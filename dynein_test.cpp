#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

#define EPSILON 0.05

int equal(float v1, float v2) {
	return fabs(v1 - v2)/fabs(v1) <= EPSILON;
}

int main() {
	
	float one, two;
	
	//dyn->set_bla( ( 90 / 180) * M_PI );
	//dyn->set_mla( ( 90 / 180) * M_PI );
	//dyn->set_mra( (-90 / 180) * M_PI );
	//dyn->set_bra( (-90 / 180) * M_PI );
	
	forces f;
	f.r_blx = 0;     f.r_bly = 0;
	f.r_mlx = 0;     f.r_mly = 0;
	f.r_tx  = 0;	 f.r_ty  = 0;
	f.r_mrx = 0;     f.r_mry = 0;
	f.r_brx = 0;     f.r_bry = 0;

	f.f_blx = 0;     f.f_bly = 0;
	f.f_mlx = 0;     f.f_mly = 0;
	f.f_tx  = 0;     f.f_ty  = 0;
	f.f_mrx = 0;     f.f_mry = 0;
	f.f_brx = 0;     f.f_bry = 0;
	
	printf("\n\nRunning dynein_struct member value check...\n\n");
	
	Dynein* dyn = new Dynein(( 90 / 180) * M_PI, ( 90 / 180) * M_PI, ( -90 / 180) * M_PI, ( -90 / 180) * M_PI, f);
	
	printf("\tIs brx defined properly..........");
	one = dyn->get_brx();
	two = -lt*sin(dyn->get_bla())*dyn->get_d_bla() + -ls*sin(dyn->get_mla())*dyn->get_d_mla() + ls*sin(dyn->get_mra())*dyn->get_d_mra() + lt*sin(dyn->get_bra()) * dyn->get_d_bra();
	if (equal(one, two))
	{
		printf("yes.\n");
	} else 
	{
		printf("no! Exiting.\n\n");
		exit(EXIT_FAILURE);
	}
	
	exit(EXIT_SUCCESS);
}













