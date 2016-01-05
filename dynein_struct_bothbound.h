#include <math.h>
#include "MersenneTwister.h"

typedef struct
{
  double nbx;   double nby;
  double nmx;   double nmy;
  double tx;    double ty;
  double fmx;   double fmy;
  double fbx;   double fby;
} bothbound_forces;

typedef struct
{
  double nba, nma, ta, fma, fba;
} bothbound_equilibrium_angles;

const bothbound_equilibrium_angles bothbound_pre_powerstroke_internal_angles = {
  0.6 * M_PI,
  0.6 * M_PI,
  0.6 * M_PI,
  1.4 * M_PI,
  0.4 * M_PI
};

/* ******************* BOTHBOUND DYNEIN CLASS DEFINITION ********************** */

class Dynein_bothbound {
public:
  Dynein(double nma_init, double fma_init, double nbx_init, double nby_init,
	 bothbound_forces* internal_test, bothbound_forces* brownian_test,
	 bothbound_equilibrium_angles* eq_angles);

  void set_nma(double d);
  void set_fma(double d);
  void set_L(double d);

  double get_nma(); // actual coordinates
  double get_fma();

  double get_nba(); // utility fns, calculated from nma, fma
  double get_fba(); 
  
  double get_nbx();
  double get_nmx();
  double get_tx();
  double get_fmx();
  double get_fbx();

  double get_nby();
  double get_nmy();
  double get_ty();
  double get_fmy();
  double get_fby();
  
  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_nma();  // bothbound
  double get_d_fma();

  double get_d_nbx();
  double get_d_nmx();
  double get_d_tx();
  double get_d_fmx();
  double get_d_fbx();

  double get_d_nby();
  double get_d_nmy();
  double get_d_ty();
  double get_d_fmy();
  double get_d_fby();

  forces get_internal();
  forces get_brownian();

  double get_binding_rate();
  double get_unbinding_rate();

  double get_PE();
  double get_KE();

  void update_velocities();
  
private:
  void update_brownian_forces();
  void update_internal_forces();

  bothbound_equilibrium_angles eq;      //Equilibrium angles

  double nma, fma; //Bothbound coordinates
  double nbx, nby;
  double L;

  double d_ln;   //Bothbound velocities
  double d_lf;

  bothbound_forces r; //Brownian forces
  bothbound_forces f; //Internal Forces

  Mode mode;
  bothbound_forces brownian_testcase;
  bothbound_forces internal_testcase;
};

