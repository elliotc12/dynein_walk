#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dynein_struct.h"

int main() {
	
	Dynein* dyn = new Dynein(bla_init, mla_init, mra_init, bra_init);
	
	printf("\n\nC Leftbound Accelerations:\n\n");
	printf("d_bla: %E\n", dyn->get_d_bla());
	printf("d_mla: %E\n", dyn->get_d_mla());
	printf("d_mra: %E\n", dyn->get_d_mra());
	printf("d_bra: %E\n", dyn->get_d_bra());
	
	free(dyn);
	dyn = NULL;
	
}


