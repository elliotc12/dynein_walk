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

void resetLog() {
  FILE* data_file = fopen("data.txt", "w");
  FILE* config_file = fopen("config.txt", "w");

  fprintf(config_file, "#gb\tgm\tgt\tdt\truntime?\tstate\n");
  fprintf(config_file, "%g\t%g\t%g\t%g\t%g\n",
          (double) gb, (double) gm, (double) gt, dt, runtime);
  fprintf(data_file,
	  "#KE\tPE\tEnergy\tt\tc1x\tc1y\tc2x\tc2y\tc3x\tc3y\tc4x\tc4y\tc5x\tc5y\tS\n");

  fclose(data_file);
  fclose(config_file);
}
