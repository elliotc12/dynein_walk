#include <math.h>
#include "MersenneTwister.h"

const double lt = 10.0;
const double ls = 10.0;

const double g = 1.0;

const double inctime = 0.1;

const double bla_0 = (108.0 / 180) * M_PI;
const double  la_0 = (108.0 / 180) * M_PI;
const double  ta_0 = (108.0 / 180) * M_PI;
const double  ra_0 = (252.0 / 180) * M_PI;

typedef enum
{
  TEST_NO_INTERNAL_FORCES,
  TEST_LEFT_INTERNAL_FORCES,
  TEST_RIGHT_INTERNAL_FORCES,
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

/* ******************************** DYNEIN CLASS DEFINITION *************************************** */

class Dynein {
public:
  Dynein(double bla_init, double mla_init, double mra_init, double bra_init,
         State s, Mode m, forces *brownian_test = 0);

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
   
  double get_f_blx();
  double get_f_mlx();
  double get_f_tx();
  double get_f_mrx();
  double get_f_brx();
  
  double get_f_bly();
  double get_f_mly();
  double get_f_ty();
  double get_f_mry();
  double get_f_bry();
  
  forces get_brownian();

  double get_PE();
  double get_KE();

  State get_state();

  void log(double t);
  void resetLog();
  void update_velocities();
  void read_init_file();
  
private:
  void update_brownian_forces();
  void update_internal_forces();

  MTRand rand;

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

  double f_blx;   double f_bly;	//Internal Forces
  double f_mlx;   double f_mly;
  double f_tx;    double f_ty;
  double f_mrx;   double f_mry;
  double f_brx;   double f_bry;

  Mode mode;
  forces *brownian_testcase;
  // internal_testcase *forces;
  State state;
};

/* *********************************** UTILITY PROTOTYPES ****************************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
