#include "dynein_struct.h"

#include <stdlib.h>
#include <fstream>

extern double bla_init;
extern double mla_init;
extern double mra_init;
extern double bra_init;

double dist(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2-x1, 2) + pow(y2-y1, 2));
}

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1) * range;
}

double square(double num) {
	return num * num;
}

double cube(double num) {
	return num * num * num;
}

double fourth(double num) {
	return num * num * num * num;
}

double fifth(double num) {
	return num * num * num * num * num;
}
