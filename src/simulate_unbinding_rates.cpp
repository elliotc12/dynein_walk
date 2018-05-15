#include <cstdio>
#include <cstring>
#include <limits>

#include <errno.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "default_parameters.h"
#include "dynein_struct.h"
#include "simulation_defaults.h"

const bool display_step_info = false;
const bool display_progress = false;

extern movie_data_struct* on_crash_old_movie_data_global_ptr;
extern movie_data_struct* on_crash_new_movie_data_global_ptr;
extern char* crash_movie_file_name_global;

const bool frozen_in_bothbound = true;

bool am_making_movie = true;
bool am_debugging_onebound = false;

long num_movie_writes = 1e7;
// bytes per movie write: 213, 2000Mb bytes max movie size

static const long MAX_FILESIZE_PERMITTED = 1<<30; // <- bit shift

void on_crash_write_movie_buffer();

struct job_msg_t {
  long long max_iteration;
  double start_time;
  int write_rate;
  char *run_msg;
  FILE *data_file;
};

void log_unbinding_probability_data(FILE* data_file, void* dyn, long long iteration, int write_rate) {
  if (iteration %  (int) (1.0/(write_rate*dt)) == 0) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    fprintf(data_file, "%13.10g %10.2g %10.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g %7.2g\n",
	    iteration*dt,
	    dyn_bb->get_near_unbinding_rate(),
	    dyn_bb->get_far_unbinding_rate(),
	    dyn_bb->get_nma(),
	    dyn_bb->get_fma(),
	    dyn_bb->get_nbx(),
	    dyn_bb->get_nmx(),
	    dyn_bb->get_tx(),
	    dyn_bb->get_fmx(),
	    dyn_bb->get_fbx(),
	    dyn_bb->get_nby(),
	    dyn_bb->get_nmy(),
	    dyn_bb->get_ty(),
	    dyn_bb->get_fmy(),
	    dyn_bb->get_fby());
    fflush(data_file);
  }
}

void unbinding_probability_callback(void* dyn, State s, void *job_msg_, data_union *job_data, long long iteration) {
  job_msg_t job_msg = *(job_msg_t *)job_msg_;

  if (iteration % 1000 == 0 and ftell(stdout) > MAX_FILESIZE_PERMITTED) {
      printf("ERROR: stdout has gotten too big.  Exciting!\n");
      exit(1);
  }

  log_unbinding_probability_data(job_msg.data_file, dyn, iteration, job_msg.write_rate);

  if (job_msg.max_iteration == 0 and iteration % (int)5e5 == 0) fflush(job_msg.data_file);
}

void set_input_variables(int argc, char** argv, char* run_name, double* runtime, int* write_rate, double* L_init) {
  char c;
  char* label;
  *run_name = 0;
  label = 0;

  static struct option long_options[] =
    {
      {"ls",                          required_argument,    0, 'a'},
      {"lt",                          required_argument,    0, 'b'},
      {"cb",                          required_argument,    0, 'c'},
      {"cm",                          required_argument,    0, 'd'},
      {"ct",                          required_argument,    0, 'e'},
      {"seed",                        required_argument,    0, 'h'},
      {"k_b",                         required_argument,    0, 'i'},
      {"k_ub",                        required_argument,    0, 'j'},
      {"runtime",                     required_argument,    0, 'l'},
      {"label",                       required_argument,    0, 'm'},
      {"dt",                          required_argument,    0, 'n'},
      {"eqb",                         required_argument,    0, 'p'},
      {"eqmpre",                      required_argument,    0, 'q'},
      {"eqmpost",                     required_argument,    0, 'r'},
      {"eqt",                         required_argument,    0, 's'},
      {"write_rate",                  required_argument,    0, 't'},
      {"c",                           required_argument,    0, 'u'},
      {"L",                           required_argument,    0, 'v'},
      {0, 0, 0, 0}
    };

  int option_index = 0;

  while ((c = getopt_long(argc, argv, "a:b:c:d:e:f:", long_options, &option_index)) != -1) {
    switch (c) {
    case 0:
      if (long_options[option_index].flag != 0) // option set a flag
	break;
      else {
	printf ("unknown option %s", long_options[option_index].name);
	if (optarg)
	  printf (" with arg %s", optarg);
	printf ("\n");
	break;
      }
    case 'a':
      Ls = strtod(optarg, NULL);
      break;
    case 'b':
      Lt = strtod(optarg, NULL);
      break;
    case 'c':
      cb = strtod(optarg, NULL);
      break;
    case 'd':
      cm = strtod(optarg, NULL);
      break;
    case 'e':
      ct = strtod(optarg, NULL);
      break;
    case 'h':
      RAND_INIT_SEED = atoi(optarg);
      break;
    case 'i':
      low_affinity_binding_rate = strtod(optarg, NULL);
      break;
    case 'j':
      low_affinity_unbinding_rate = strtod(optarg, NULL);
      break;
    case 'l':
      *runtime = strtod(optarg, NULL);
      break;
    case 'm':
      label = optarg;
      break;
    case 'n':
      dt = strtod(optarg, NULL);
      break;
    case 'p':
      onebound_post_powerstroke_internal_angles.bba = strtod(optarg, NULL) * M_PI / 180.0;
      bothbound_pre_powerstroke_internal_angles.nba = strtod(optarg, NULL) * M_PI / 180.0;
      bothbound_pre_powerstroke_internal_angles.fba = strtod(optarg, NULL) * M_PI / 180.0;
      break;
    case 'q': // pre-powerstroke
      onebound_post_powerstroke_internal_angles.uma = strtod(optarg, NULL) * M_PI / 180.0;
      break;
    case 'r': // post-powerstroke
      onebound_post_powerstroke_internal_angles.bma = strtod(optarg, NULL) * M_PI / 180.0;
      bothbound_pre_powerstroke_internal_angles.nma = strtod(optarg, NULL) * M_PI / 180.0;
      bothbound_pre_powerstroke_internal_angles.fma = strtod(optarg, NULL) * M_PI / 180.0;
      break;
    case 's':
      onebound_post_powerstroke_internal_angles.ta = strtod(optarg, NULL) * M_PI / 180.0;
      bothbound_pre_powerstroke_internal_angles.ta = strtod(optarg, NULL) * M_PI / 180.0;
      break;
    case 't':
      *write_rate = strtod(optarg, NULL);
      break;
    case 'u':
      exponential_unbinding_angle_constant = strtod(optarg, NULL);
      break;
    case 'v':
      *L_init = strtod(optarg, NULL);
      break;
    case '?':
      printf("Some other unknown getopt error.\n");
      exit(EXIT_FAILURE);
    default:
      printf("Default case in getopt: uh-oh!\n");
      exit(EXIT_FAILURE);
    }
  }

  sprintf(run_name, "%s__L-%g,s-%g", label, *L_init, RAND_INIT_SEED);

  if (optind != argc) {
    printf("Improper usage!");
    exit(EXIT_FAILURE);
  }
}

void sig_handler_print_movie_buffer(int signum) {
  if (am_making_movie) {
    printf("Received sigusr1, writing to data file and turning off printing.\n");
    on_crash_write_movie_buffer();
    am_making_movie = false;
    am_only_writing_on_crash = false;
  }
  else {
    printf("Turning on printing mode.\n");
    am_making_movie = true;
    am_only_writing_on_crash = true;
  }
}

int main(int argc, char** argv) {
  fflush(stdout);
  setvbuf(stdout, 0, _IONBF, 0);

  char* run_name = new char[200];
  double runtime = 0;
  int write_rate = 0;
  double L_init = 0;
  set_input_variables(argc, argv, run_name, &runtime, &write_rate, &L_init);


  char *data_fname = new char[200];
  sprintf(data_fname, "data/unbinding_probability/%s.txt", run_name);

  job_msg_t job_msg;
  job_msg.max_iteration = 0;
  job_msg.start_time = clock();
  job_msg.write_rate = write_rate;
  job_msg.run_msg = run_name;
  job_msg.data_file = fopen(data_fname, "w");

  fprintf(job_msg.data_file, "# command line:");
  for (int i=0; i<argc; i++) {
    fprintf(job_msg.data_file, " %s", argv[i]);
  }

  printf("\n\n\n*********%s*********\n", run_name);
  fprintf(job_msg.data_file, "\n\n\n\n#********%s********\n", run_name);
  fprintf(job_msg.data_file, "#%12s, %9s, %9s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s, %6s\n",
	  "t", "near_p_ub", "far_p_ub", "nma", "fma", "nbx", "nmx", "tx", "fmx", "fbx", "nby", "nmy", "ty", "fmy", "fby");

  if (errno) {
    perror("Error opening stepping data or movie file.\n");
    exit(errno);
  }

  double init_position[] = {120.0*M_PI / 180.0, 120.0*M_PI / 180.0, 0, 0, L_init}; // nma_init, fma_init, nbx, nby, L

  printf("Initial conditions: nma %g fma %g L %g\n", init_position[0], init_position[1], init_position[4]);

  simulate(runtime, RAND_INIT_SEED, BOTHBOUND, init_position, unbinding_probability_callback, &job_msg, NULL);

  fclose(job_msg.data_file);

  return EXIT_SUCCESS;
}
