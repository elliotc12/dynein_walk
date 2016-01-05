#include <math.h>
#include "MersenneTwister.h"

typedef struct
{
  double bbx;   double bby;
  double bmx;   double bmy;
  double tx;    double ty;
  double umx;   double umy;
  double ubx;   double uby;
} onebound_forces;

typedef struct
{
  double bba, bma, ta, uma;
} onebound_equilibrium_angles;

const onebound_equilibrium_angles onebound_post_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.0 * M_PI,
  0.6 * M_PI
};

/* ******************** ONEBOUND DYNEIN CLASS DEFINITION ********************** */

class Dynein_onebound {
public:
  Dynein(double bba_init, double bma_init, double fma_init, double fba_init,
	       double bbx_init, double bby_init, State s, forces *internal_test,
	       forces *brownian_test, equilibrium_angles* eq_angles);

  /** Onebound functions **/
  void set_bba(double d);
  void set_bma(double d);
  void set_uma(double d);
  void set_uba(double d);

  void set_bbx(double d);
  void set_bby(double d);
  
  double get_bba();
  double get_bma();
  double get_uma();
  double get_uba();

  double get_bbx();
  double get_bmx();
  double get_tx();
  double get_umx();
  double get_ubx();

  double get_bby();
  double get_bmy();
  double get_ty();
  double get_umy();
  double get_uby();

  void set_state(State s);
  
  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_bba();
  double get_d_bma();
  double get_d_uma();
  double get_d_uba();

  double get_d_bbx();
  double get_d_bmx();
  double get_d_tx();
  double get_d_umx();
  double get_d_ubx();

  double get_d_bby();
  double get_d_bmy();
  double get_d_ty();
  double get_d_umy();
  double get_d_uby();

  forces get_internal();
  forces get_brownian();

  double get_binding_rate();
  double get_unbinding_rate();

  double get_PE();
  double get_KE();

  State get_state();
  
  void update_velocities();
  
private:
  void update_brownian_forces();
  void update_internal_forces();

  onebound_equilibrium_angles eq;      //Equilibrium angles

  double bba;    //Onebound coordinates
  double bma;
  double uma;
  double uba;
  
  double bbx, bby;
  
  double d_bba;   //Onebound angular velocities
  double d_bma;
  double d_uma;
  double d_uba;

  onebound_forces r; //Brownian forces
  onebound_forces f; //Internal Forces

  Mode mode;
  onebound_forces *brownian_testcase;
  onebound_forces *internal_testcase;
  State state;
};
