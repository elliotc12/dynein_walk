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

void write_onebound_data_callback(void* dyn, State s, void** job_msg, data_union *job_data, long long iteration) {
  if (iteration % pe_calculation_skip_iterations != 0) return;
  else {
    assert(s == NEARBOUND);
    long long max_iteration = *((long long**) job_msg)[0];
    double start_time = *((double**) job_msg)[1];
    char* run_msg = ((char**) job_msg)[2];

    FILE* data_file = ((FILE**) job_msg)[3];

    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;

    onebound_data_generate_struct data;
    data.time = iteration*dt;
    data.bba_PE = dyn_ob->PE_bba;
    data.bma_PE = dyn_ob->PE_bma;
    data.ta_PE  = dyn_ob->PE_ta;
    data.uma_PE = dyn_ob->PE_uma;
    data.bba = dyn_ob->get_bba();
    data.bma = dyn_ob->get_bma() + M_PI - dyn_ob->get_bba();
    data.ta  = dyn_ob->get_uma() - dyn_ob->get_bma();
    data.uma = dyn_ob->get_uma() + M_PI - dyn_ob->get_uba();

    fwrite(&data, sizeof(onebound_data_generate_struct), 1, data_file);

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

  char data_fname[200];

  strcpy(data_fname, "data/onebound_data_");
  strcpy(config_fname, "data/ob_config_");

  strcat(data_fname, f_appended_name);
  strcat(config_fname, f_appended_name);

  strcat(data_fname, ".bin");
  strcat(config_fname, ".txt");

  prepare_data_file(NULL, data_fname);
  write_config_file(config_fname, 0, "Initial state: onebound\nInitial conformation: equilibrium\n");

  void* job_msg[4];
  job_msg[0] = (double*) &iterations;

  double current_time = clock();
  job_msg[1] = &current_time;

  char run_msg[512];
  sprintf(run_msg, "seed = %d", (int) RAND_INIT_SEED);
  job_msg[2] = run_msg;

  FILE* data_file = fopen(data_fname, "a");

  struct stat data_stat;
  if (stat(data_fname, &data_stat) == -1) {
    perror("Error using stat on data_file.bin");
    exit(1);
  }

  if (setvbuf(data_file, NULL, _IOFBF, data_stat.st_blksize) == -1) {
    perror("Error using setvbuf");
    exit(1);
  }

  job_msg[3] =  data_file;

  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  double init_position[] = {eq.bba, eq.bma, eq.ta, eq.uma, 0, 0};

  simulate(iterations*dt, RAND_INIT_SEED, NEARBOUND, init_position,
	   write_onebound_data_callback, job_msg, NULL);

  fclose(data_file);
  return EXIT_SUCCESS;
}
