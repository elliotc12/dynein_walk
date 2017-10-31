#include <cassert>
#include <fenv.h>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <errno.h>

#include "dynein_struct.h"
#include "simulation_defaults.h"
#include "plotting_defaults.h"

extern double bla_init;
extern double mla_init;
extern double mra_init;
extern double bra_init;

char* crash_movie_file_name_global;

movie_data_struct* on_crash_old_movie_data_global_ptr; // ugly, fix sometime
movie_data_struct* on_crash_new_movie_data_global_ptr;

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

void prepare_data_file(const char* custom_str, char* fname) {
  if (custom_str == NULL) {
    FILE* data_file = fopen(fname, "w"); // clear the file
    if (errno) {
      perror("Error preparing data file");
      printf("data file: %s\n", fname);
      exit(errno);
    }
    fclose(data_file);
  }
  else {
    FILE* data_file = fopen(fname, "w");
    if (errno) {
      perror("Error preparing data file");
      printf("data file: %s\n", fname);
      exit(errno);
    }
    fprintf(data_file, "%s\n", custom_str);
    fclose(data_file);
  }
}

void append_data_to_file(double* data1, double* data2, int len, FILE* file) {
    for (int i = 0; i < len; i++) {
      assert(data1[i] == data1[i]);
      // assert(data2[i] == data2[i]);
      fprintf(file, "%+.12e     %+.5e\n", data1[i], data2[i]);
    }
}

void write_config_file(char* fname, int omit_flags, const char* custom_str) {
  char text_buf[2048];
  char buf[100];
  text_buf[0] = 0;
  if (custom_str != NULL) strcat(text_buf, custom_str);
  sprintf(buf, "Lt: %g\n", Lt);
  strcat(text_buf, buf);
  sprintf(buf, "Ls: %g\n", Ls);
  strcat(text_buf, buf);
  sprintf(buf, "T: %g\n", T);
  if ((omit_flags & CONFIG_OMIT_T) == 0) strcat(text_buf, buf);
  sprintf(buf, "ct: %.2e\n", ct);
  if ((omit_flags & CONFIG_OMIT_C) == 0) strcat(text_buf, buf);
  sprintf(buf, "cm: %.2e\n", cm);
  if ((omit_flags & CONFIG_OMIT_C) == 0) strcat(text_buf, buf);
  sprintf(buf, "cb: %.2e\n", cb);
  if ((omit_flags & CONFIG_OMIT_C) == 0) strcat(text_buf, buf);
  sprintf(buf, "dt: %g\n", dt);
  strcat(text_buf, buf);
  if ((omit_flags & CONFIG_INCLUDE_SKIPINFO) != 0) {
    if (custom_generate_averaging_width != 0) {
      sprintf(buf, "l. avg width: %d\n", custom_generate_averaging_width);
      strcat(text_buf, buf);
    }
    else {
      sprintf(buf, "l. avg width: between pts\n");
      strcat(text_buf, buf);
    }
    sprintf(buf, "l. avg pe points: %d\n", num_generate_pe_datapoints);
    strcat(text_buf, buf);
    sprintf(buf, "l. avg angle points: %d\n", num_generate_angle_datapoints);
    strcat(text_buf, buf);
  }

  FILE* data_file = fopen(fname, "w");
  fputs(text_buf, data_file);
  fclose(data_file);
}

void write_movie_config(char* movie_config_fname, double runtime) {
  FILE* config_file = fopen(movie_config_fname, "w");
  fprintf(config_file,
	  "#gb\t"
	  "gm\t"
	  "gt\t"
	  "dt\t"
	  "runtime?\t"
	  "state\t"
	  "kbT\n");
  fprintf(config_file, "%g\t%g\t%g\t%g\t%g\t%g\n",
          (double) gb, (double) gm, (double) gt, dt, runtime, kb*T);
  fclose(config_file);
}

void on_crash_write_movie_buffer() {
  movie_data_struct* old_data = on_crash_old_movie_data_global_ptr;
  movie_data_struct* new_data = on_crash_new_movie_data_global_ptr;
  FILE* data_file = fopen(crash_movie_file_name_global, "w");

  fprintf(data_file, "#State\ttime\tPE_1\tPE_2\tPE_3\tPE_4\tPE_5\t"
	  "x1\ty1\tx2\ty2\tx3\ty3\tx4\ty4\tx5\ty5\t"
	  "fx1\tfy1\tfx2\tfy2\tfx3\tfy3\tfx4\tfy4\tfx5\tfy5\n");

  int i = 0;
  while (old_data[i].time > 0 and i < MOVIE_BUFFER_SIZE) {
    fprintf(data_file,
	    "%d\t"
	    "%.10g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t"
	    "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
	    "\n",
	    old_data[i].state,
	    old_data[i].time,
	    old_data[i].PE_1, old_data[i].PE_2, old_data[i].PE_3, old_data[i].PE_4, old_data[i].PE_5,
	    old_data[i].x_1, old_data[i].y_1, old_data[i].x_2, old_data[i].y_2, old_data[i].x_3,
	    old_data[i].y_3, old_data[i].x_4, old_data[i].y_4, old_data[i].x_5, old_data[i].y_5,
	    old_data[i].fx_1, old_data[i].fy_1, old_data[i].fx_2, old_data[i].fy_2, old_data[i].fx_3,
	    old_data[i].fy_3, old_data[i].fx_4, old_data[i].fy_4, old_data[i].fx_5, old_data[i].fy_5);
    i++;
  }
  int j = 0;
  while (new_data[j].time > 0) {
    fprintf(data_file,
	    "%d\t"
	    "%.10g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t%g\t%g\t"
	    "%g\t%g\t"
	    "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g"
	    "\n",
	    new_data[j].state,
	    new_data[j].time,
	    new_data[j].PE_1, new_data[j].PE_2, new_data[j].PE_3, new_data[j].PE_4, new_data[j].PE_5,
	    new_data[j].x_1,  new_data[j].y_1,  new_data[j].x_2,  new_data[j].y_2,  new_data[j].x_3,
	    new_data[j].y_3,  new_data[j].x_4,  new_data[j].y_4,  new_data[j].x_5,  new_data[j].y_5,
	    new_data[j].fx_1, new_data[j].fy_1, new_data[j].fx_2, new_data[j].fy_2, new_data[j].fx_3,
	    new_data[j].fy_3, new_data[j].fx_4, new_data[j].fy_4, new_data[j].fx_5, new_data[j].fy_5);
    j++;
  }
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

DynArr::DynArr(int init_len) {
  len = init_len;
  current = 0;
  data = (double*)malloc(sizeof(double)*len);
  if (data == NULL) {
    perror("error allocating memory for DynArr");
    exit(errno);
  }
}

DynArr::~DynArr() {
  free(data);
}

void DynArr::append(double d) {
  if (current == len) {
    len *= 2;
    data = (double*)realloc(data, sizeof(double)*len);
    if (data == NULL) {
      perror("error reallocating memory for DynArr");
      exit(errno);
    }
  }
  data[current] = d;
  current += 1;
}

double* DynArr::get_data() {
  return data;
}
