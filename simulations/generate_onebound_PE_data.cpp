#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <limits>

#include "../default_parameters.h"
#include "../dynein_struct.h"

void append_pe_data_to_file(double time, double pe, FILE* file) {
  pe_data data;
  data.time = time;
  data.pe = pe;
  fwrite(&data, sizeof(pe_data), 1, file);
}

void write_onebound_PEs_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  int skip_num = *((int**) job_msg)[7];
  if (iteration % skip_num != 0)
    return;
  else {
    assert(s == NEARBOUND);
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];

    FILE* bb_file = ((FILE**) job_msg)[3];
    FILE* bm_file = ((FILE**) job_msg)[4];
    FILE*  t_file = ((FILE**) job_msg)[5];
    FILE* um_file = ((FILE**) job_msg)[6];

    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;

    double time = iteration*dt;

    append_pe_data_to_file(time, dyn_ob->PE_bba, bb_file);
    append_pe_data_to_file(time, dyn_ob->PE_bma, bm_file);
    append_pe_data_to_file(time, dyn_ob->PE_ta ,  t_file);
    append_pe_data_to_file(time, dyn_ob->PE_uma, um_file);

    if (iteration % (skip_num*10) == 0) {
      printf("PE calculation progress for %s: %lld / %lld, %g%%                \r",
    	     run_msg, iteration, max_iteration, ((double) iteration) / max_iteration * 100);
      fflush(NULL);
    }

    if (iteration == (max_iteration-1) - ((max_iteration-1) % skip_num)) { // final iteration
      printf("Finished generating PE data for %s, which took %g seconds               \n", run_msg,
	     ((double) clock() - start_time) / CLOCKS_PER_SEC);
    }
  }
}

int main(int argc, char** argv) {
  MICROTUBULE_BINDING_DISTANCE = -std::numeric_limits<double>::infinity();
  MICROTUBULE_REPULSION_FORCE = 0.0;

  long long iterations = 1e11;
  int skip_num = 1e4;

  if (argc != 2) {
    printf("Error, TITLE variable must have underscores, not spaces.\n");
    exit(1);
  }

  char* f_appended_name = argv[1];
  char bb_fname[200];
  char bm_fname[200];
  char t_fname[200];
  char um_fname[200];
  char config_fname[200];

  strcpy(bb_fname, "data/bba_pe_vs_time_");
  strcpy(bm_fname, "data/bma_pe_vs_time_");
  strcpy(t_fname, "data/ta_pe_vs_time_");
  strcpy(um_fname, "data/uma_pe_vs_time_");
  strcpy(config_fname, "data/config_pe_vs_time_");

  strcat(bb_fname, f_appended_name);
  strcat(bm_fname, f_appended_name);
  strcat(t_fname, f_appended_name);
  strcat(um_fname, f_appended_name);
  strcat(config_fname, f_appended_name);

  strcat(bb_fname, ".bin");
  strcat(bm_fname, ".bin");
  strcat(t_fname, ".bin");
  strcat(um_fname, ".bin");
  strcat(config_fname, ".txt");

  prepare_data_file(NULL, bb_fname);
  prepare_data_file(NULL, bm_fname);
  prepare_data_file(NULL,  t_fname);
  prepare_data_file(NULL, um_fname);
  write_config_file(config_fname, 0, NULL);

  void* job_msg[8];
  job_msg[0] = (double*) &iterations;

  double current_time = clock();
  job_msg[1] = &current_time;

  const char* run_msg_base = "pe gathering (";
  char run_msg[512];
  strcpy(run_msg, run_msg_base);
  char seedbuf[50];
  sprintf(seedbuf, "seed = %d)", (int) RAND_INIT_SEED);
  strcat(run_msg, seedbuf);
  job_msg[2] = run_msg;

  FILE* bb_file = fopen(bb_fname, "a");
  FILE* bm_file = fopen(bm_fname, "a");
  FILE*  t_file = fopen( t_fname, "a");
  FILE* um_file = fopen(um_fname, "a");

  job_msg[3] = bb_file;
  job_msg[4] = bm_file;
  job_msg[5] = t_file;
  job_msg[6] = um_file;
  job_msg[7] = &skip_num;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, init_position,
	   write_onebound_PEs_callback, job_msg, NULL);

  fclose(bb_file);
  fclose(bm_file);
  fclose( t_file);
  fclose(um_file);
  return EXIT_SUCCESS;
}
