#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

#define EPSILON 0.05

int equal(float v1, float v2) {
	if (v1 == 0) { return fabs(v2) < EPSILON; }
	else return fabs(v1 - v2)/fabs(v1) <= EPSILON;
}

int main() {
	
	float one, two;
	
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
	
	//Dynein* dyn = new Dynein(( 90 / 180)*M_PI, ( 90 / 180)*M_PI, ( -90 / 180)*M_PI, ( -90 / 180)*M_PI, f);
	Dynein* dyn = new Dynein((108.0 / 180) * M_PI, (36.0 / 180) * M_PI, -(144.0 / 180) * M_PI, -(72.0 / 180) * M_PI, f);
	
	printf("\tIs d_blx defined properly..........");
	one = dyn->get_d_blx();
	two = 0;
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_bly defined properly..........");
	one = dyn->get_d_bly();
	two = 0;
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_mlx defined properly..........");
	one = dyn->get_d_mlx();
	two = -ls*sin(dyn->get_bla())*dyn->get_d_bla();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_mly defined properly..........");
	one = dyn->get_d_mly();
	two = ls*cos(dyn->get_bla())*dyn->get_d_bla();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs  d_tx defined properly..........");
	one = dyn->get_d_tx();
	two = -ls*sin(dyn->get_bla())*dyn->get_d_bla() + -lt*sin(dyn->get_mla())*dyn->get_d_mla();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs  d_ty defined properly..........");
	one = dyn->get_d_ty();
	two = ls*cos(dyn->get_bla())*dyn->get_d_bla() + lt*cos(dyn->get_mla())*dyn->get_d_mla();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_mrx defined properly..........");
	one = dyn->get_d_mrx();
	two = -ls*sin(dyn->get_bla())*dyn->get_d_bla() + -lt*sin(dyn->get_mla())*dyn->get_d_mla() + lt*sin(dyn->get_mra())*dyn->get_d_mra();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_mry defined properly..........");
	one = dyn->get_d_mry();
	two = ls*cos(dyn->get_bla())*dyn->get_d_bla() + lt*cos(dyn->get_mla())*dyn->get_d_mla() + -lt*cos(dyn->get_mra())*dyn->get_d_mra();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_brx defined properly..........");
	one = dyn->get_d_brx();
	two = -ls*sin(dyn->get_bla())*dyn->get_d_bla() + -lt*sin(dyn->get_mla())*dyn->get_d_mla() + lt*sin(dyn->get_mra())*dyn->get_d_mra() + ls*sin(dyn->get_bra())*dyn->get_d_bra();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	printf("\tIs d_bry defined properly..........");
	one = dyn->get_d_bry();
	two = ls*cos(dyn->get_bla())*dyn->get_d_bla() + lt*cos(dyn->get_mla())*dyn->get_d_mla() + -lt*cos(dyn->get_mra())*dyn->get_d_mra() + -ls*cos(dyn->get_bra())*dyn->get_d_bra();
	if (equal(one, two))
	{
		printf("yes, struct says %f and calculation says %f.\n", one, two);
	} else 
	{
		printf("no! struct says %f and calculation says %f. Exiting.\n\n", one, two);
		exit(EXIT_FAILURE);
	}
	
	exit(EXIT_SUCCESS);
}













