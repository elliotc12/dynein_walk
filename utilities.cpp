#include <cassert>
#include <fenv.h>
#include <fstream>
#include <stdlib.h>

#include "dynein_struct.h"

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

double get_average(double* data, int len) {
  assert(len != 0);
  double sum = 0;
  for (int i = 0; i < len; i++) {
    sum += data[i];
  }
  return sum / len;
}

double get_variance(double* data, int len) {
  assert(len != 0);
  double sum = 0;
  double ave = get_average(data, len);
  for (int i = 0; i < len; i++) {
    sum += (data[i] - ave)*(data[i] - ave);
  }
  return sum / len;
}

void FPE_signal_handler(int signum) {
  fexcept_t flag;
  fegetexceptflag(&flag, FE_ALL_EXCEPT);
  printf("FE exception occured.\n");
  printf("FE_INVALID | flag: %d\n", FE_INVALID & flag);
  printf("FE_DIVBYZERO | flag: %d\n", FE_DIVBYZERO & flag);
  printf("FE_INEXACT | flag: %d\n", FE_INEXACT & flag);
  printf("FE_UNDERFLOW | flag: %d\n", FE_UNDERFLOW & flag);
  printf("FE_OVERFLOW | flag: %d\n", FE_OVERFLOW & flag);
  exit(EXIT_FAILURE);
}

#ifdef __APPLE__
void feenableexcept(int x) { printf("fake feenableexcept for mac.\n"); }
#endif
