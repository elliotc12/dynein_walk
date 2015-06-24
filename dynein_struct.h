#include <math.h>

const double lt = 10.0;
const double ls = 10.0;

const double g = 1;

const double kt  = 1.0; // Higher spring constant -> less deviation from equilibrium
const double kml = 1.0;
const double kmr = 1.0;
const double kbl = 1.0;
const double kbr = 1.0;

const double mb = 1.0; // Higher mass constant -> less movement
const double mm = 1.0;
const double mt = 1.0;

const double ba = (108.0 / 180) * M_PI;
const double ma = (108.0 / 180) * M_PI;
const double ta = (108.0 / 180) * M_PI;

const double bla_init = (108.0 / 180) * M_PI;
const double mla_init = (36.0 / 180) * M_PI;
const double mra_init = (36.0 / 180) * M_PI;
const double bra_init = (108.0 / 180) * M_PI;

const double inctime = 0.1;
const double runtime = 100.0;

extern double fblx, fbly;		// Brownian variables
extern double fmlx, fmly;
extern double ftx, fty;
extern double fmrx, fmry;
extern double fbrx, fbry;

typedef enum {
	LEFTBOUND,
	RIGHTBOUND,
	BOTHBOUND
} states;


/* ******************************** DYNEIN CLASS DEFINITION *************************************** */

class Dynein {
public:
  Dynein(double bla_init, double mla_init, double mra_init, double bra_init);	
	
  void set_bla(double d);
  void set_mla(double d);
  void set_mra(double d);
  void set_bra(double d);

  void set_blx(double d);
  void set_bly(double d);

  void set_state(states s);

  void set_d_bla(double d);
  void set_d_mla(double d);
  void set_d_mra(double d);
  void set_d_bra(double d);

  double get_bla();
  double get_mla();
  double get_mra();
  double get_bra();

  double get_blx() {
    return blx;
  }
  double get_mlx();
  double get_tx();
  double get_mrx();
  double get_brx();

  double get_bly();
  double get_mly();
  double get_ty();
  double get_mry();
  double get_bry();
  
  // The following are dynamical properties that only exist in an
  // ephemeral per-timestep way:
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

  double get_d_bla();
  double get_d_mla();
  double get_d_mra();
  double get_d_bra();
   
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
  
  double get_r_blx();
  double get_r_mlx();
  double get_r_tx();
  double get_r_mrx();
  double get_r_brx();
  
  double get_r_bly();
  double get_r_mly();
  double get_r_ty();
  double get_r_mry();
  double get_r_bry();

  double get_PE();
  double get_KE();

  states get_state();

  void log(double t);
  void update_protein();
  
private:
  double bla;
  double mla;
  double mra;
  double bra;

  double blx, bly;
  double brx;     double bry;

  double mlx;     double mly;
  double tx;      double ty;
  double mrx;     double mry;
  
  double d_bla;   //Angular Velocities
  double d_mla;
  double d_mra;
  double d_bra;

  double d_blx;   double d_bly;
  double d_mlx;   double d_mly;
  double d_tx;    double d_ty;
  double d_mrx;   double d_mry;
  double d_brx;   double d_bry;

  double r_blx;   double r_bly; //Brownian forces
  double r_mlx;   double r_mly;
  double r_tx;    double r_ty;
  double r_mrx;   double r_mry; 					
  double r_brx;   double r_bry;					

  double f_blx;   double f_bly;	//Internal Forces
  double f_mlx;   double f_mly;
  double f_tx;    double f_ty;
  double f_mrx;   double f_mry;
  double f_brx;   double f_bry;

  states state;
};

/* *********************************** UTILITY PROTOTYPES ****************************************** */
double randAngle(double range);
double dist(double d, double h, double i, double j);
void   resetLog(Dynein* dyn);
double square(double num);
double cube(double num);
double fourth(double num);
double fifth(double num);
