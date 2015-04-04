#include "dynein_struct.h"

#include <stdlib.h>
#include <fstream>


double dist(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1) * range;
}

void resetLog(Dynein* dyn) {
	FILE* data_file = fopen("data.txt", "w");
	FILE* config_file = fopen("config.txt", "w");
	
	fprintf(config_file, "#inctime\truntime\tstate\n%+.3f\t%+.3f\t%d\n", inctime, runtime, (int) dyn->get_state());
	fprintf(data_file,
		"#KE\t\t\t\tPE\t\t\t\tEnergy\t\tt\t\tblX\t\t\tblY\t\t\tmlX\t\t\tmlY\t\t\ttX\t\t\ttY\t\t\tmrX\t\t\tmrY\t\t\tbrX\t\t\tbrY\t\t\tS\n");
	
	fclose(data_file);
	fclose(config_file);
}

double square(double num) {
	return num * num;
}

double cube(double num) {
	return num * num * num;
}

double fifth(double num) {
	return num * num * num * num * num;
}
