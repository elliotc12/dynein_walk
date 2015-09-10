#include <math.h>
#include "MersenneTwister.h"

const double kb = 0.1;
const double T = 1.0;

const double lt = 10.0;
const double ls = 10.0;

const double gm = 2.0;  // motor domain gamma
const double gb = 1.0;  // binding domain gamma
const double gt = 3.0;  // tail domain gamma

const double inctime = 0.1;

typedef enum
{
  PRE_POWERSTROKE,
  POST_POWERSTROKE
} Mode;

typedef enum
{
  LEFTBOUND,
  RIGHTBOUND,
  BOTHBOUND
} State;

typedef struct
{
  double blx;   double bly;
  double mlx;   double mly;
  double tx;    double ty;
  double mrx;   double mry;
  double brx;   double bry;
} forces;

typedef struct
{
  double bla, la, ta, ra;
} equilibrium_angles;

const equilibrium_angles pre_powerstroke_leftbound_internal_angles = {
  (108.0 / 180) * M_PI,
  (108.0 / 180) * M_PI,
  (108.0 / 180) * M_PI,
  (252.0 / 180) * M_PI
};

const equilibrium_angles pre_powerstroke_rightbound_internal_angles = {
  (72.0 / 180) * M_PI,
  (252.0 / 180) * M_PI,
  (252.0 / 180) * M_PI,
  (108.0 / 180) * M_PI
};

/* ******************************** DYNEIN CLASS DEFINITION *************************************** */

class Dynein {
public:
  Dynein(double bla_init, double mla_init, double mra_init, double bra_init,
         State s, forces* internal_test, forces* brownian_test, equilibrium_angles* eq_angles);

  void set_bla(double d);
  void set_mla(double d);
  void set_mra(double d);
  void set_bra(double d);

  void set_blx(double d);
  void set_bly(double d);

  void set_state(State s);

  double get_bla();
  double get_mla();
  double get_mra();
  double get_bra();

  double get_blx();
  double get_mlx();
  double get_tx();
  double get_mrx();
  double get_brx();
                  ;
  double get_bly();
  double get_mly();
  double get_ty();
  double get_mry();
  double get_bry();
  
  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:

  double get_d_bla();
  double get_d_mla();
  double get_d_mra();
  double get_d_bra();
  
  double get_d_blx();
  double get_d_mlx();
  double get_d_tx();
  double get_d_mrx();
  double get_d_brx();
  
  double get_d_bly();
  double get_d_mly();
  double get_d_ty();
  double get_d_mry();
  double get_d_bry();
   
  forces get_internal();
  forces get_brownian();

  double get_PE();
  double get_KE();

  State get_state();

  void log(double t);
  void resetLog();
  void update_velocities();
  
private:
  void update_brownian_forces();
  void update_internal_forces();

  MTRand rand;

  equilibrium_angles eq;

  double bla;
  double mla;
  double mra;
  double bra;
  
  double blx, bly;
  
  double d_bla;   //Angular Velocities
  double d_mla;
  double d_mra;
  double d_bra;

  forces r; //Brownian forces
  forces f; //Internal Forces

  Mode mode;
  forces *brownian_testcase;
  forces *internal_testcase;
  State state;
};

/* *********************************** UTILITY PROTOTYPES ****************************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
