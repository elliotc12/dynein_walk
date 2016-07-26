#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include <sys/stat.h>
#include <unistd.h>

#include "../default_parameters.h"
#include "../dynein_struct.h"
#include "simulation_defaults.h"

void append_pe_data_to_file(double time, double pe, FILE* file) {
  data_2d data;
  data.time = time;
  data.d = pe;
  fwrite(&data, sizeof(data_2d), 1, file);
}

void write_onebound_data_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  if (iteration % pe_calculation_skip_iterations != 0)
    return;
  else {
    assert(s == NEARBOUND);
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];

    FILE* bb_pe_file = ((FILE**) job_msg)[3];
    FILE* bm_pe_file = ((FILE**) job_msg)[4];
    FILE*  t_pe_file = ((FILE**) job_msg)[5];
    FILE* um_pe_file = ((FILE**) job_msg)[6];

    FILE* bb_angle_file = ((FILE**) job_msg)[7];
    FILE* bm_angle_file = ((FILE**) job_msg)[8];
    FILE*  t_angle_file = ((FILE**) job_msg)[9];
    FILE* um_angle_file = ((FILE**) job_msg)[10];

    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;

    double time = iteration*dt;

    double bba = dyn_ob->get_bba();
    double bma = dyn_ob->get_bma() + M_PI - dyn_ob->get_bba();
    double ta  = dyn_ob->get_uma() - dyn_ob->get_bma();
    double uma = dyn_ob->get_uma() + M_PI - dyn_ob->get_uba();

    append_pe_data_to_file(time, dyn_ob->PE_bba, bb_pe_file);
    append_pe_data_to_file(time, dyn_ob->PE_bma, bm_pe_file);
    append_pe_data_to_file(time, dyn_ob->PE_ta ,  t_pe_file);
    append_pe_data_to_file(time, dyn_ob->PE_uma, um_pe_file);

    append_pe_data_to_file(time, bba, bb_angle_file);
    append_pe_data_to_file(time, bma, bm_angle_file);
    append_pe_data_to_file(time, ta ,  t_angle_file);
    append_pe_data_to_file(time, uma, um_angle_file);

    if (iteration % (pe_calculation_skip_iterations*10) == 0) {
      printf("PE calculation progress (%s): %lld / %lld, %g%%                \r", run_msg,
    	     iteration, max_iteration, ((double) iteration) / max_iteration * 100);
      fflush(NULL);
    }

    if (iteration == (max_iteration-1) - ((max_iteration-1) % pe_calculation_skip_iterations)) {
      printf("Finished generating PE data (%s), process took %g seconds               \n", run_msg,
	     ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
}

int main(int argc, char** argv) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;

  long long iterations = 1e12;
  T = 500;

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char* f_appended_name = argv[1];
  char config_fname[200];

  char bb_pe_fname[200];
  char bm_pe_fname[200];
  char t_pe_fname[200];
  char um_pe_fname[200];

  char bb_angle_fname[200];
  char bm_angle_fname[200];
  char t_angle_fname[200];
  char um_angle_fname[200];

  strcpy(bb_pe_fname, "data/bba_pe_data_");
  strcpy(bm_pe_fname, "data/bma_pe_data_");
  strcpy(t_pe_fname, "data/ta_pe_data_");
  strcpy(um_pe_fname, "data/uma_pe_data_");
  strcpy(bb_angle_fname, "data/bba_angle_data_");
  strcpy(bm_angle_fname, "data/bma_angle_data_");
  strcpy(t_angle_fname, "data/ta_angle_data_");
  strcpy(um_angle_fname, "data/uma_angle_data_");
  strcpy(config_fname, "data/config_");

  strcat(bb_pe_fname, f_appended_name);
  strcat(bm_pe_fname, f_appended_name);
  strcat(t_pe_fname, f_appended_name);
  strcat(um_pe_fname, f_appended_name);
  strcat(bb_angle_fname, f_appended_name);
  strcat(bm_angle_fname, f_appended_name);
  strcat(t_angle_fname, f_appended_name);
  strcat(um_angle_fname, f_appended_name);
  strcat(config_fname, f_appended_name);

  strcat(bb_pe_fname, ".bin");
  strcat(bm_pe_fname, ".bin");
  strcat(t_pe_fname, ".bin");
  strcat(um_pe_fname, ".bin");
  strcat(bb_angle_fname, ".bin");
  strcat(bm_angle_fname, ".bin");
  strcat(t_angle_fname, ".bin");
  strcat(um_angle_fname, ".bin");
  strcat(config_fname, ".txt");

  prepare_data_file(NULL, bb_pe_fname);
  prepare_data_file(NULL, bm_pe_fname);
  prepare_data_file(NULL,  t_pe_fname);
  prepare_data_file(NULL, um_pe_fname);
  prepare_data_file(NULL, bb_angle_fname);
  prepare_data_file(NULL, bm_angle_fname);
  prepare_data_file(NULL,  t_angle_fname);
  prepare_data_file(NULL, um_angle_fname);
  write_config_file(config_fname, 0, "Initial state: onebound\nInitial conformation: equilibrium\n");

  void* job_msg[11];
  job_msg[0] = (double*) &iterations;

  double current_time = clock();
  job_msg[1] = &current_time;

  const char* run_msg_base = "";
  char run_msg[512];
  strcpy(run_msg, run_msg_base);
  char seedbuf[50];
  sprintf(seedbuf, "seed = %d", (int) RAND_INIT_SEED);
  strcat(run_msg, seedbuf);
  job_msg[2] = run_msg;

  FILE* bb_pe_file = fopen(bb_pe_fname, "a");
  FILE* bm_pe_file = fopen(bm_pe_fname, "a");
  FILE*  t_pe_file = fopen( t_pe_fname, "a");
  FILE* um_pe_file = fopen(um_pe_fname, "a");
  FILE* bb_angle_file = fopen(bb_angle_fname, "a");
  FILE* bm_angle_file = fopen(bm_angle_fname, "a");
  FILE*  t_angle_file = fopen( t_angle_fname, "a");
  FILE* um_angle_file = fopen(um_angle_fname, "a");

  struct stat bb_stat;
  if (stat(bb_pe_fname, &bb_stat) == -1) {
    perror("Error using stat on bb_file");
    exit(1);
  }

  if (setvbuf(bb_pe_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf(bm_pe_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf( t_pe_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf(um_pe_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf(bb_angle_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf(bm_angle_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf( t_angle_file, NULL, _IOFBF, bb_stat.st_blksize) == -1 or
      setvbuf(um_angle_file, NULL, _IOFBF, bb_stat.st_blksize) == -1) {
    perror("Error using setvbuf");
    exit(1);
  }

  job_msg[3] =  bb_pe_file;
  job_msg[4] =  bm_pe_file;
  job_msg[5] =  t_pe_file;
  job_msg[6] =  um_pe_file;
  job_msg[7] =  bb_angle_file;
  job_msg[8] =  bm_angle_file;
  job_msg[9] =  t_angle_file;
  job_msg[10] = um_angle_file;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, init_position,
	   write_onebound_data_callback, job_msg, NULL);

  fclose(bb_pe_file);
  fclose(bm_pe_file);
  fclose( t_pe_file);
  fclose(um_pe_file);
  fclose(bb_angle_file);
  fclose(bm_angle_file);
  fclose( t_angle_file);
  fclose(um_angle_file);
  return EXIT_SUCCESS;
}
