#include "../default_parameters.h"
#include "../dynein_struct.h"

void plot_log(void* dyn, State s, void* job_msg, void* job_data, int iteration) {
  FILE* data_file = (FILE*) job_msg;
  if (s == NEARBOUND or s == FARBOUND) {
    Dynein_onebound* dyn_ob = (Dynein_onebound*) dyn;
    onebound_forces f = dyn_ob->get_internal();
    fprintf(data_file,
	    "%d\t"
	    "%.2g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
	    "\n",
	    dyn_ob->get_state(),
	    iteration*dt,
	    dyn_ob->PE_bba, dyn_ob->PE_bma, dyn_ob->PE_ta, dyn_ob->PE_uma, 0.0,
	    dyn_ob->get_bbx(), dyn_ob->get_bby(),
	    dyn_ob->get_bmx(), dyn_ob->get_bmy(),
	    dyn_ob->get_tx(), dyn_ob->get_ty(),
	    dyn_ob->get_umx(), dyn_ob->get_umy(),
	    dyn_ob->get_ubx(), dyn_ob->get_uby(),
	    f.bbx, f.bby, f.bmx, f.bmy, f.tx, f.ty, f.umx, f.umy, f.ubx, f.uby);
  }
  else if (s == BOTHBOUND) {
    Dynein_bothbound* dyn_bb = (Dynein_bothbound*) dyn;
    bothbound_forces f = dyn_bb->get_internal();
    fprintf(data_file,
	    "%d\t"
	    "%.2g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
	    "\n",
	    BOTHBOUND,
	    iteration*dt,
	    dyn_bb->PE_nba, dyn_bb->PE_nma, dyn_bb->PE_ta, dyn_bb->PE_fma, dyn_bb->PE_fba,
	    dyn_bb->get_nbx(), dyn_bb->get_nby(),
	    dyn_bb->get_nmx(), dyn_bb->get_nmy(),
	    dyn_bb->get_tx(), dyn_bb->get_ty(),
	    dyn_bb->get_fmx(), dyn_bb->get_fmy(),
	    dyn_bb->get_fbx(), dyn_bb->get_fby(),
	    f.nbx, f.nby, f.nmx, f.nmy, f.tx, f.ty,f.fmx, f.fmy, f.fbx, f.fby);
  } 
  else { // must be unbound
    fprintf(data_file,
	    "%d\t"
	    "%.2g\t"
	    "%.2g\t%.2g\t%.2g\t%.2g\t%.2g\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t"
	    "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
	    "\n",
	    UNBOUND,
	    iteration*dt,
	    0.0, 0.0, 0.0, 0.0, 0.0,
	    0.0, 0.0,
	    0.0, 0.0,
	    0.0, 0.0,
	    0.0, 0.0,
	    0.0, 0.0,
	    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
}

void reset_logs(FILE* data_file, FILE* config_file, double runtime) {
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
  fprintf(data_file,
  	  "#State\t"
  	  "%10s\t"
  	  "PE.b1\tPE.m1\tPE.t\tPE.m2\tPE.b2\t"
  	  "b1x\tb1y\t"
  	  "m1x\tm1y\t"
  	  "tx\tty\t"
  	  "m2x\tm2y\t"
  	  "b2x\tb2y\t"
  	  "f.b1x\tf.b1y\tf.m1x\tf.m1y\tf.tx\tf.ty\tf.m2x\tf.m2y\tf.b2x\tf.b2y"
  	  "\n", "t");
}

int main(int argvc, char **argv) {
  if (argvc != 6) {
    printf("Error. Usage: ./walk runtime bla_init mla_init mra_init bra_init.\n");
    return EXIT_FAILURE;
  }

  T = 100;
  double runtime = strtod(argv[1], NULL) * dt;
  onebound_equilibrium_angles eq = onebound_post_powerstroke_internal_angles;
  
  double bba_init = strtod(argv[2], NULL) * M_PI + eq.bba;
  double bma_init = strtod(argv[3], NULL) * M_PI + bba_init + eq.bma - M_PI;
  double uma_init = strtod(argv[4], NULL) * M_PI + bma_init + eq.ta;
  double uba_init = strtod(argv[5], NULL) * M_PI + uma_init + M_PI - eq.uma;

  printf("%s %s %s %s %s %s\n", argv[0], argv[1], argv[2], argv[3], argv[4], argv[5]);
  printf("bba: %g, bma: %g, uma: %g, uba: %g\n", bba_init, bma_init, uma_init, uba_init);
  
  double init_position[] = {bba_init, bma_init, uma_init, uba_init, 0, 0};

  FILE* data_file = fopen("data.txt", "w+");
  FILE* config_file = fopen("config.txt", "w+");
  reset_logs(data_file, config_file, runtime);

  void* job_msg = (void*) data_file;
  
  simulate(runtime, RAND_INIT_SEED, NEARBOUND, init_position, plot_log, job_msg, NULL);
  
  fclose(data_file);
  fclose(config_file);
  
  return EXIT_SUCCESS;
}
