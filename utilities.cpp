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

void resetLogs(FILE* data_file, FILE* config_file) {
  fprintf(config_file, "#gb\tgm\tgt\tdt\truntime?\tstate\n");
  fprintf(config_file, "%g\t%g\t%g\t%g\t%g\n",
          (double) gb, (double) gm, (double) gt, dt, runtime);
  fprintf(data_file,
	  "#KE\tPE\tEnergy\tt\tc1x\tc1y\tc2x\tc2y\tc3x\tc3y\tc4x\tc4y\tc5x\tc5y\tS\n");
}

void log_run(FILE* run_file, double runtime, double run_length,
			      double distance_traveled, int steps) {

  float ave_step_dist = distance_traveled / steps;
  float ave_step_time = runtime / steps;

  printf("\n\n***********Run data**********\n");
  printf("Run length: %g nm\n", run_length);
  printf("Distance traveled: %g nm\n", distance_traveled);
  printf("Steps: %d\n", steps);
  printf("Average step length: %g nm\n", ave_step_dist);
  printf("Average step time: %g s\n\n\n", ave_step_time);
  fprintf(run_file, "Run length \tDistance traveled \tSteps \tAve step length \tAve step time\n");
  fprintf(run_file, "%f\t%f\t%d\t%f\t%g\n",
  	  run_length, distance_traveled, steps, ave_step_dist, ave_step_time);
  fclose(run_file);
}
