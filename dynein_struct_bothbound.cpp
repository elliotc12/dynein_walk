#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct.h"

inline double sqr(double x) {
  return x*x;
}

/* ********************* BOTHBOUND DYNEIN FUNCTIONS ****************************** */

Dynein_bothbound::Dynein_bothbound(double nma_init, double fma_init, double nbx_init,
                                   double nby_init, double L_init,
                                   bothbound_forces* internal_test,
                                   bothbound_forces* brownian_test,
                                   bothbound_equilibrium_angles* eq_angles) {
  nbx = nbx_init;
  nby = nby_init;

  nma = nma_init;
  fma = fma_init;

  L = L_init;

  internal_testcase = internal_test;
  brownian_testcase = brownian_test;

  if (eq_angles) {
    eq = *eq_angles; // use test angles
  } else {
    eq = bothbound_pre_powerstroke_internal_angles; // use experimental angles
  }

  update_velocities();
}

Dynein_bothbound::Dynein_bothbound(Dynein_onebound* old_dynein, MTRand* mtrand) { // out of old dyn
  rand = mtrand;

  if (old_dynein->get_state() == State::NEARBOUND) {
    nbx = old_dynein->get_bbx();
    nby = 0;

    nma = M_PI + old_dynein->get_bma() - old_dynein->get_bba();
    fma = M_PI + old_dynein->get_uma() - old_dynein->get_uba();

    L = old_dynein->get_ubx() - old_dynein->get_bbx();

  } else {
    nbx = old_dynein->get_ubx();
    nby = 0;

    nma = M_PI + old_dynein->get_uma() - old_dynein->get_uba();
    fma = M_PI + old_dynein->get_bma() - old_dynein->get_bba();

    L = old_dynein->get_bbx() - old_dynein->get_ubx();
  }

  update_velocities();
}

void Dynein_bothbound::update_brownian_forces() {
  if (brownian_testcase) {
    r = *brownian_testcase; // just copy over forces!
  } else {
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.nbx, &r.nby);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.nmx, &r.nmy);
    rand->gauss2(sqrt(2*kb*T/(gt*dt)), &r.tx, &r.ty);
    rand->gauss2(sqrt(2*kb*T/(gm*dt)), &r.fmx, &r.fmy);
    rand->gauss2(sqrt(2*kb*T/(gb*dt)), &r.fbx, &r.fby);
  }
}

void Dynein_bothbound::update_internal_forces() {
  if (internal_testcase) {
    f = *internal_testcase;
  } else {
    f.nbx = 0;     f.nby = 0;     // Initialize forces to zero
    f.nmx = 0;     f.nmy = 0;
    f.tx  = 0;     f.ty  = 0;
    f.fmx = 0;     f.fmy = 0;
    f.fbx = 0;     f.fby = 0;

    double T, f1, f2, f1x, f1y, f2x, f2y;

    T = cb*(get_nba() - eq.nba);
    f2 = T/ls;
    f2x = f2 * sin(get_nba());
    f2y = f2 * -cos(get_nba());
    f.nmx += f2x;
    f.nmy += f2y;
    f.nbx += -f2x; // Equal and opposite forces!  :)
    f.nby += -f2y; // Equal and opposite forces!  :)

    T = cm*(nma - eq.nma);
    f1 = T/ls;
    f2 = T/lt;
    f1x = f1 * sin(get_nba());
    f1y = f1 * -cos(get_nba());
    f2x = f2 * sin(nma - (M_PI - get_nba()));
    f2y = f2 * -cos(nma - (M_PI - get_nba()));
    f.nbx += f1x;
    f.nby += f1y;
    f.tx  += f2x;
    f.ty  += f2y;
    f.nmx += -(f1x + f2x);
    f.nmy += -(f1y + f2y);

    T = ct*(fma + get_fba() - nma - get_nba() - eq.ta);
    f1 = T / lt;
    f2 = T / lt;
    f1x = f1 * sin(nma + get_nba() - M_PI);
    f1y = f1 * -cos(nma + get_nba() - M_PI);
    f2x = f2 * -sin(fma + get_fba() - M_PI);  // not sure if these angles are right?
    f2y = f2 * cos(fma + get_fba() - M_PI);
    f.nmx += f1x;
    f.nmy += f1y;
    f.fmx += f2x;
    f.fmy += f2y;
    f.tx  += -(f1x + f2x);
    f.ty  += -(f1y + f2y);

    T = cm*(fma - eq.fma);
    f1 = T / lt;
    f2 = T / ls;
    f1x = f1 * sin(fma + get_fba() - M_PI);
    f1y = f1 * -cos(fma + get_fba() - M_PI);
    f2x = f2 * sin(get_fba());
    f2y = f2 * -cos(get_fba());
    f.tx  += f1x;
    f.ty  += f1y;
    f.fbx += f2x;
    f.fby += f2y;
    f.fmx += -(f1x + f2x);
    f.fmy += -(f1y + f2y);

    T = cb*(get_fba() - eq.fba);
    f1 = T / ls;
    f1x = f1 * sin(get_fba());
    f1y = f1 * -cos(get_fba());
    f.fmx += f1x;
    f.fmy += f1y;
    f.fbx += -f1x;
    f.fby += -f1y;

    if (get_nby() < 0) f.nby += MICROTUBULE_REPULSION_FORCE * fabs(get_nby());
    if (get_nmy() < 0) f.nmy += MICROTUBULE_REPULSION_FORCE * fabs(get_nmy());
    if (get_ty()  < 0) f.ty  += MICROTUBULE_REPULSION_FORCE * fabs(get_ty());
    if (get_fmy() < 0) f.fmy += MICROTUBULE_REPULSION_FORCE * fabs(get_fmy());
    if (get_fby() < 0) f.fby += MICROTUBULE_REPULSION_FORCE * fabs(get_fby());
  }
}

void Dynein_bothbound::update_velocities() {
  update_internal_forces();
  update_brownian_forces();

  double pm_n = 1; // plus or minus
  double pm_f = 1;
  if (nma > M_PI) pm_n = -1;
  if (fma > M_PI) pm_f = -1;

  const double Ln = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(nma));
  const double Lf = sqrt(sqr(Ls) + sqr(Lt) - 2*Ls*Lt*cos(fma));

  double Axn = (pow(L, 2)*pow(lt, 2) - pow(L, 2)*pow(ls, 2)) / (2*(L*pow(Ln, 3))) 
  + (Lf*pow(lt, 2) - Lf*pow(ls, 2)) / (2*(L*pow(Ln, 2))) 
  + (pow(Lf, 2)*pow(ls, 2) - pow(Lf, 2)*pow(lt, 2)) / (2*(L*pow(Ln, 3))) 
  - ls*(pm_n*((pow(L, 2) / (2*pow(Ln, 3)) 
  + (pow(Lf, 4) - 2*(pow(L, 2)*pow(Lf, 2))) / (2*(pow(L, 2)*pow(Ln, 3))) 
  - Ln / (2*pow(L, 2)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Ln, 2)*pow(ls, 2))) 
         + (2*pow(lt, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))) 
      / (2*sqrt(1 - pow(L, 2) / (4*pow(Ln, 2)) 
  + (2*(pow(L, 2)*pow(Lf, 2)) - pow(Lf, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
  + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2)))) 
  - ls*(pm_n*(((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2))) 
  / (2*(pow(Ln, 3)*pow(ls, 2))) - Ln / (2*pow(ls, 2)))
 *sqrt(1 - pow(L, 2) / (4*pow(Ln, 2)) 
         + (2*(pow(L, 2)*pow(Lf, 2)) - pow(Lf, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
         + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
      / (2*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(Ln, 2)*pow(ls, 2))) 
  + (2*pow(lt, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1));

  double Axf = ls*(pm_n*(((Lf*pow(L, 2) - pow(Lf, 3)) / (pow(L, 2)*pow(Ln, 2)) + Lf / pow(L, 2))
                         *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4))
                               / (4*(pow(Ln, 2)*pow(ls, 2)))
                               + (2*pow(lt, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)))
    / (-2*sqrt(1 - pow(L, 2) / (4*pow(Ln, 2))
               + (2*(pow(L, 2)*pow(Lf, 2)) - pow(Lf, 4)) / (4*(pow(L, 2)*pow(Ln, 2)))
               + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))));
  
  double Bxn = Ln/L;
  
  double Bxf = -Lf/L;

  double Cxn = - pm_f * ls *
    sqrt(1 - (pow(ls,4) + pow(lt,4) - 2*sqr(ls)*sqr(lt)
              )/(4*pow(ls,2)*pow(Lf,2)) -
         (pow(Lf,2) + 2*pow(ls,2) - 2*pow(lt,2)
          )/(4*pow(ls,2)))
    /(2*sqrt(1 - pow(L,2)/(4*pow(Lf,2))
             - (Ln^4 - 2*pow(L,2)*pow(Ln,2)
                )/(4*pow(L,2)*pow(Lf,2))
             - (pow(Lf,2) + 2*pow(L,2) - 2*pow(Ln,2)
                )/(4*pow(L,2))))
    (-(Ln^3 - pow(L,2)*Ln
       )/(pow(L,2)*pow(Lf,2)) + Ln/pow(L,2)); // check if this is right

  double Cxf = (pow(L, 2)*pow(lt, 2) - pow(L, 2)*pow(ls, 2)) / (2*(L*pow(Lf, 3))) 
  + (Ln*pow(lt, 2) - Ln*pow(ls, 2)) / (2*(L*pow(Lf, 2))) 
  + (pow(Ln, 2)*pow(ls, 2) - pow(Ln, 2)*pow(lt, 2)) / (2*(L*pow(Lf, 3))) 
  - ls*(pm_f*((pow(L, 2) / (2*pow(Lf, 3)) 
  + (pow(Ln, 4) - 2*(pow(L, 2)*pow(Ln, 2))) / (2*(pow(L, 2)*pow(Lf, 3))) 
  - Lf / (2*pow(L, 2)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Lf, 2)*pow(ls, 2))) 
         + (2*pow(lt, 2) - pow(Lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))) 
      / (2*sqrt(1 - pow(L, 2) / (4*pow(Lf, 2)) 
  + (2*(pow(L, 2)*pow(Ln, 2)) - pow(Ln, 4)) / (4*(pow(L, 2)*pow(Lf, 2))) 
  + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2)))) 
  - ls*(pm_f*(((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2))) 
  / (2*(pow(Lf, 3)*pow(ls, 2))) - Lf / (2*pow(ls, 2)))
 *sqrt(1 - pow(L, 2) / (4*pow(Lf, 2)) 
         + (2*(pow(L, 2)*pow(Ln, 2)) - pow(Ln, 4)) / (4*(pow(L, 2)*pow(Lf, 2))) 
         + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
      / (2*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(Lf, 2)*pow(ls, 2))) 
  + (2*pow(lt, 2) - pow(Lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1));
  
  double Ayn = ((pow(lt, 2) - pow(ls, 2) - pow(Ln, 2)) / (2*pow(Ln, 2)) + 1)
 *sqrt(1 - pow(L, 2) / (4*pow(Ln, 2)) 
         + (2*(pow(L, 2)*pow(Lf, 2)) - pow(Lf, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
         + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))) 
  + (pow(L, 2) / (2*pow(Ln, 3)) 
       + (pow(Lf, 4) - 2*(pow(L, 2)*pow(Lf, 2))) / (2*(pow(L, 2)*pow(Ln, 3))) 
       - Ln / (2*pow(L, 2)))*(pow(ls, 2) + pow(Ln, 2) - pow(lt, 2)) 
      / (4*(Ln*sqrt(1 - pow(L, 2) / (4*pow(Ln, 2)) 
    + (2*(pow(L, 2)*pow(Lf, 2)) - pow(Lf, 4)) / (4*pow(L, 2)*pow(Ln,2))
  + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm_n*(ls*(1 / L + (pow(Lf, 2) - pow(L, 2) - pow(Ln, 2)) / (2*pow(L*Ln, 2)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Ln, 2)*pow(ls, 2))) 
         + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + ls*(pm_n*((pow(L, 2) + pow(Ln, 2) - pow(Lf, 2))
 *((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2))) 
     / (2*(pow(Ln, 3)*pow(ls, 2))) - Ln / (2*pow(ls, 2))))) 
      / (4*(L*Ln
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Ln, 2)*pow(ls, 2))) 
         + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)));
  
  double Ayf = ((Lf*pow(L, 2) - pow(Lf, 3)) / (pow(L, 2)*pow(Ln, 2)) + Lf / pow(L, 2))
 *(pow(ls, 2) + pow(Ln, 2) - pow(lt, 2)) 
  / (4*(Ln*sqrt(1 - pow(L, 2) / (4*pow(Ln, 2)) 
  + (2*(pow(L, 2)*pow(Lf, 2)) - pow(Lf, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
  + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm_n*(ls*(-(Lf / (L*Ln)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Ln, 2)*pow(ls, 2))) 
         + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + Lf*(pm_n*(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2))) 
      / (4*(L*(Ln*(ls*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(Ln, 2)*pow(ls, 2))) 
  + (2*pow(Lf, 2) - pow(Ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)))));
  
  double Byn = (2*Ln - Ln*(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / pow(L, 2)) 
  / (2*sqrt(pow(Ln, 2) 
  - pow(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2), 2) / (4*pow(L, 2))));
  
  double Byf = Lf*(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) 
  / (2*(pow(L, 2)
 *sqrt(pow(Ln, 2) - pow(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2), 2) / (4*pow(L, 2)))));
  
  double Cyn = ((Ln*pow(L, 2) - pow(Ln, 3)) / (pow(L, 2)*pow(Ln, 2)) + Ln / pow(L, 2))
 *(pow(ls, 2) + pow(Lf, 2) - pow(lt, 2)) 
  / (4*(Lf*sqrt(1 - pow(L, 2) / (4*pow(Lf, 2)) 
  + (2*(pow(L, 2)*pow(Ln, 2)) - pow(Ln, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
  + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm_f*(ls*(-(Ln / (L*Lf)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Lf, 2)*pow(ls, 2))) 
         + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + Ln*(pm_f*(pow(L, 2) + pow(Lf, 2) - pow(Ln, 2))) 
      / (4*(L*(Lf*(ls*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(Lf, 2)*pow(ls, 2))) 
  + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)))));
  
  double Cyf = ((pow(lt, 2) - pow(ls, 2) - pow(Lf, 2)) / (2*pow(Lf, 2)) + 1)
 *sqrt(1 - pow(L, 2) / (4*pow(Lf, 2)) 
         + (2*(pow(L, 2)*pow(Ln, 2)) - pow(Ln, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
         + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))) 
  + (pow(L, 2) / (2*pow(Lf, 3)) 
       + (pow(Ln, 4) - 2*(pow(L, 2)*pow(Ln, 2))) / (2*(pow(L, 2)*pow(Lf, 3))) 
       - Lf / (2*pow(L, 2)))*(pow(ls, 2) + pow(Lf, 2) - pow(lt, 2)) 
      / (4*(Lf*sqrt(1 - pow(L, 2) / (4*pow(Lf, 2)) 
  + (2*(pow(L, 2)*pow(Ln, 2)) - pow(Ln, 4)) / (4*(pow(L, 2)*pow(Ln, 2))) 
  + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm_f*(ls*(1 / L + (pow(Ln, 2) - pow(L, 2) - pow(Lf, 2)) / (2*(L*pow(Lf, 2))))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(Lf, 2)*pow(ls, 2))) 
         + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + ls*(pm_f*((pow(L, 2) + pow(Lf, 2) - pow(Ln, 2))
 *((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2)))
     / (2*(pow(Lf, 3)*pow(ls, 2))) - Lf / (2*pow(ls, 2))))) 
      / (4*(L*(Lf*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(Lf, 2)*pow(ls, 2))) 
  + (2*pow(Ln, 2) - pow(Lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))));
  
  d_Ln = (Cyf*f.fmx - Cxf*f.fmy + Cyf*g*r.fmx - Cxf*g*r.fmy)/((Cxn*Cyf - Cxf*Cyn)*g) - (((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy)) +
    (Cxn*Cyf - Cxf*Cyn)*(-((Byn*f.fmx)/g) + (Cxn*f.ty)/g - Byn*r.fmx + Cxn*r.ty))*(get_nmx() - get_tx()))/g -
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn) *(-((Bxn*f.fmx)/g) + (Cxn*f.tx)/g - Bxn*r.fmx +
    Cxn*r.tx) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy))*(get_nmy() -
    get_ty()))/g)*(Cyf*get_fmx() - Cyf*get_tx() - Cxf*get_fmy() + Cxf*get_ty()))/ ((Cxn*Cyf - Cxf*Cyn)*g*(-((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(get_nmy() - get_ty())*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(get_fmx() - get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/g) -
    (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(get_nmx() - get_tx())*((Cxn*Cyf - Cxf*Cyn)* ((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) -
    (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn* (-get_fmy() + get_ty()))/g)))/g)) - (((-(Cyf*get_fbx()) +
    Cyf*get_fmx() + Cxf*get_fby() - Cxf*get_fmy())/((Cxn*Cyf - Cxf*Cyn)*g) - (((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() -
    get_tx())*((Byn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g -
    (Cxn*(get_fby() - get_fmy()))/g)))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Bxn*(Cxn* Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g -
    (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g))*(get_nmy() - get_ty()))/g)*(Cyf*get_fmx() -
    Cyf*get_tx() - Cxf*get_fmy() + Cxf*get_ty()))/((Cxn*Cyf - Cxf*Cyn)*g*(-((Cxn* (Cxn*Cyf - Cxf*Cyn)*(get_nmy() -
    get_ty())*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(get_fmx() - get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/ g) - (-(Bxn*Cxf) +
    Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() -
    get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g)))* ((-((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Ayn*f.fmx)/g) + (Cxn*f.nmy)/g -
    Ayn*r.fmx + Cxn*r.nmy) - (-(Ayn*Cxf) + Ayf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx +
    Cxn*r.fmy))*(get_nbx() - get_nmx()))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*((Cxn*f.nmx)/g -
    (Axn*f.fmx)/g + Cxn*r.nmx - Axn*r.fmx) - (-(Axn*Cxf) + Bxf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g -
    Cyn*r.fmx + Cxn*r.fmy))* (get_nby() - get_nmy()))/g)*(get_nmy() - get_ty()))/g) - (-((-(Byn*Cxf) +
    Byf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy)) + (Cxn*Cyf -
    Cxf*Cyn)*(-((Byn*f.fmx)/g) + (Cxn*f.ty)/g - Byn*r.fmx + Cxn*r.ty))* (-(((Cxn*Cxn)*((Cxn*Cyf -
    Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() - get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf -
    Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/(g*g)))* (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy()
    - get_ty())*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(get_fmx() - get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/g) - (-(Bxn*Cxf) +
    Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() -
    get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) - ((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy)) +
    (Cxn*Cyf - Cxf*Cyn)*(-((Byn*f.fmx)/g) + (Cxn*f.ty)/g - Byn*r.fmx + Cxn*r.ty))* (get_nmx() - get_tx()))/g -
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Bxn*f.fmx)/g) + (Cxn*f.tx)/g - Bxn*r.fmx +
    Cxn*r.tx) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy))*(get_nmy() -
    get_ty()))/g)* (-((-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() -
    get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() +
    get_ty()))/(g*g))* ((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g))) - (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() -
    get_ty())*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nbx() - get_nmx())*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() + get_tx()))/g -
    (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/ g) + (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(get_nby() - get_nmy())*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() + get_tx()))/g - (-(Axn*Cxf) +
    Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g))/g)))/ ((-((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nbx() - get_nmx())*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g -
    (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g)))/g) + (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(get_nby() - get_nmy())*((Axn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Axn*Cxf) +
    Bxf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g)))/g)*(get_nmy() - get_ty()))/g) - ((Byn*(Cxn*Cyf -
    Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g))*
    (-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() - get_nmy()))/(g*g)) +
    ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/(g*g)))*
    (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() - get_ty())*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(get_fmx() - get_tx()))/g) + (Bxn*(-get_fmx()
    + get_tx()))/g) - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) +
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() - get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() -
    get_ty()))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) -
    ((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() - get_tx())*((Byn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g)))/g - (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*((Bxn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g
    - (Cxn*(get_fby() - get_fmy()))/g))*(get_nmy() - get_ty()))/g)* (-((-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf -
    Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() - get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf -
    Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/(g*g))* ((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy()
    - get_ty()))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g))) - (Cxn*(Cxn*Cyf
    - Cxf*Cyn)*(get_nmy() - get_ty())*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nbx() - get_nmx())*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() +
    get_tx()))/g - (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/ g) +
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nby() - get_nmy())*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() + get_tx()))/g - (-(Axn*Cxf) +
    Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g))/g)); // from Mathematica

  
  d_Lf = (-(Cyn*f.fmx) + Cxn*f.fmy - Cyn*g*r.fmx + Cxn*g*r.fmy)/((Cxn*Cyf - Cxf*Cyn)*g) -
           (((Cxn*(Cxn*Cyf - Cxf*Cyn)* (-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g -
           Cyn*r.fmx + Cxn*r.fmy)) + (Cxn*Cyf - Cxf*Cyn)*(-((Byn*f.fmx)/g) + (Cxn*f.ty)/g - Byn*r.fmx +
           Cxn*r.ty))*(get_nmx() - get_tx()))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Bxn*f.fmx)/g)
           + (Cxn*f.tx)/g - Bxn*r.fmx + Cxn*r.tx) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*f.fmx)/g) +
           (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy))* (get_nmy() - get_ty()))/g)*(-(Cyn*get_fmx()) + Cyn*get_tx() + Cxn*get_fmy() -
           Cxn*get_ty()))/ ((Cxn*Cyf - Cxf*Cyn)*g*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() - get_ty())*((Cxn*Cyf -
           Cxf*Cyn)*(-((Cxn*(get_fmx() - get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/g) - (-(Bxn*Cxf) +
           Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx()
           - get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) +
           Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g)) - (((Cyn*get_fbx() - Cyn*get_fmx() -
           Cxn*get_fby() + Cxn*get_fmy())/((Cxn*Cyf - Cxf*Cyn)*g) - (((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() - get_tx())*
           ((Byn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g -
           (Cxn*(get_fby() - get_fmy()))/g)))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Bxn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() -
           get_fmx()))/g - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g))*(get_nmy() -
           get_ty()))/g)* (-(Cyn*get_fmx()) + Cyn*get_tx() + Cxn*get_fmy() - Cxn*get_ty()))/ ((Cxn*Cyf -
           Cxf*Cyn)*g*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() - get_ty())*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(get_fmx() -
           get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/g) - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy()
           + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() - get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() +
           get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() +
           get_ty()))/g)))/g)))* ((-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf -
           Cxf*Cyn)*(-((Ayn*f.fmx)/g) + (Cxn*f.nmy)/g - Ayn*r.fmx + Cxn*r.nmy) - (-(Ayn*Cxf) +
           Ayf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy))*(get_nbx() - get_nmx()))/g) +
           (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*((Cxn*f.nmx)/g - (Axn*f.fmx)/g + Cxn*r.nmx -
           Axn*r.fmx) - (-(Axn*Cxf) + Bxf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx +
           Cxn*r.fmy))* (get_nby() - get_nmy()))/g)*(get_nmy() - get_ty()))/g) - (-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*f.fmx)/g) +
           (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy)) + (Cxn*Cyf - Cxf*Cyn)*(-((Byn*f.fmx)/g) + (Cxn*f.ty)/g
           - Byn*r.fmx + Cxn*r.ty))* (-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() +
           get_tx())*(get_nby() - get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() -
           get_nmx())*(-get_nmy() + get_ty()))/(g*g)))* (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() - get_ty())*((Cxn*Cyf -
           Cxf*Cyn)*(-((Cxn*(get_fmx() - get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/g) - (-(Bxn*Cxf) +
           Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx()
           - get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) +
           Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g) - ((Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*f.fmx)/g) + (Cxn*f.fmy)/g - Cyn*r.fmx +
           Cxn*r.fmy)) + (Cxn*Cyf - Cxf*Cyn)*(-((Byn*f.fmx)/g) + (Cxn*f.ty)/g - Byn*r.fmx + Cxn*r.ty))*
           (get_nmx() - get_tx()))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Bxn*f.fmx)/g) +
           (Cxn*f.tx)/g - Bxn*r.fmx + Cxn*r.tx) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*f.fmx)/g) +
           (Cxn*f.fmy)/g - Cyn*r.fmx + Cxn*r.fmy))*(get_nmy() - get_ty()))/g)* (-((-(((Cxn*Cxn)*((Cxn*Cyf -
           Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() - get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf -
           Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() + get_ty()))/(g*g))* ((Cxn*Cyf -
           Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx()
           + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g))) - (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() - get_ty())*(-((Cxn*(Cxn*Cyf
           - Cxf*Cyn)*(get_nbx() - get_nmx())*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() + get_tx()))/g - (-(Ayn*Cxf) +
           Ayf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/ g) + (Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(get_nby() - get_nmy())*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() + get_tx()))/g - (-(Axn*Cxf) +
           Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g)))/g))/g)))/ ((-((Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nbx() - get_nmx())*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g
           - (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g)))/g) + (Cxn*(Cxn*Cyf
           - Cxf*Cyn)*(get_nby() - get_nmy())*((Axn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Axn*Cxf) +
           Bxf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g)))/g)*(get_nmy() - get_ty()))/g) - ((Byn*(Cxn*Cyf
           - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() -
           get_fmy()))/g))* (-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() -
           get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() +
           get_ty()))/(g*g)))* (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmy() - get_ty())*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(get_fmx() -
           get_tx()))/g) + (Bxn*(-get_fmx() + get_tx()))/g) - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy()
           + get_ty()))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() - get_tx())*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() +
           get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() +
           get_ty()))/g)))/g) - ((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nmx() - get_tx())*((Byn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() -
           get_fmx()))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g)))/g -
           (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Bxn*(Cxn*Cyf - Cxf*Cyn)*(get_fbx() - get_fmx()))/g - (-(Bxn*Cxf) +
           Bxf*Cxn)*((Cyn*(get_fbx() - get_fmx()))/g - (Cxn*(get_fby() - get_fmy()))/g))*(get_nmy() - get_ty()))/g)*
           (-((-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-get_nmx() + get_tx())*(get_nby() -
           get_nmy()))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(get_nbx() - get_nmx())*(-get_nmy() +
           get_ty()))/(g*g))* ((Cxn*Cyf - Cxf*Cyn)*((Byn*(-get_fmx() + get_tx()))/g - (Cxn*(get_fmy() - get_ty()))/g) -
           (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() + get_ty()))/g))) - (Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(get_nmy() - get_ty())*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nbx() - get_nmx())*((Ayn*(Cxn*Cyf -
           Cxf*Cyn)*(-get_fmx() + get_tx()))/g - (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy() +
           get_ty()))/g)))/ g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(get_nby() - get_nmy())*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-get_fmx() +
           get_tx()))/g - (-(Axn*Cxf) + Bxf*Cxn)*((Cyn*(-get_fmx() + get_tx()))/g - (Cxn*(-get_fmy()                         + get_ty()))/g)))/g))/g));
  // from Mathematica
}

double Dynein_bothbound::get_near_unbinding_rate() {
  if (f.nby + r.nby >= BOTHBOUND_UNBINDING_FORCE) {
    return 1.0;
  } else return 0.0;
}

double Dynein_bothbound::get_far_unbinding_rate() {
  if (f.fby + r.fby >= BOTHBOUND_UNBINDING_FORCE) {
    return 1.0;
  } else return 0.0;
}

/*** Set positions, velocities and forces ***/
void Dynein_bothbound::set_nma(double d) {
    nma = d;
}

void Dynein_bothbound::set_fma(double d) {
    fma = d;
}

void Dynein_bothbound::set_L(double d) {
    L = d;
}

/*** Angular Velocities ***/

double Dynein_bothbound::get_d_bba() {
    return d_bba;
}

double Dynein_bothbound::get_d_bma() {
    return d_bma;
}

double Dynein_bothbound::get_d_fma() {
    return d_fma;
}

double Dynein_bothbound::get_d_fba() {
    return d_fba;
}

double Dynein_bothbound::get_d_nma() {
  int pm = (nma > M_PI) ? -1 : 1; // sign of d_nma depends on value of nma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - Ln*Ln) / (2*lt*ls),2))
    * (Ln / (lt*ls)) * d_Ln;
}

double Dynein_bothbound::get_d_fma() {
  int pm = (fma > M_PI) ? -1 : 1; // sign of d_fma depends on value of fma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - Lf*Lf) / (2*lt*ls),2))
    * (Lf / (lt*ls)) * d_Lf;
}

/*** Get coordinates ***/
double Dynein_bothbound::get_nbx() {
  return nbx;
}

double Dynein_bothbound::get_nmx() {
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / (2*(L*Ln)))
         + pm*acos((pow(ls, 2) + pow(Ln, 2) - pow(lt, 2)) / (2*(Ln*ls))));
}

double Dynein_bothbound::get_tx() {
  return (pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / (2*L);
}

double Dynein_bothbound::get_fmx() {
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(Lf, 2) - pow(Ln, 2)) / (2*(L*Lf)))
		+ pm*acos((pow(ls, 2) + pow(Lf, 2) - pow(lt, 2)) / (2*(Lf*ls))));;
}

double Dynein_bothbound::get_fbx() {
  return L;
}

double Dynein_bothbound::get_nby() {
  return 0;
}

double Dynein_bothbound::get_nmy() {
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(Ln, 2) - pow(Lf, 2)) / (2*(L*Ln)))
		+ pm*acos((pow(ls, 2) + pow(Ln, 2) - pow(lt, 2)) / (2*(Ln*ls))));
}

double Dynein_bothbound::get_ty() {
  return Ln*sqrt(1 - pow(pow(L, 2) + pow(Ln, 2) - pow(Lf, 2), 2)
		 / (4*(pow(L, 2)*pow(Ln, 2))));
}

double Dynein_bothbound::get_fmy() {
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(Lf, 2) - pow(Ln, 2)) / (2*(L*Lf)))
		+ pm*acos((pow(ls, 2) + pow(Lf, 2) - pow(lt, 2)) / (2*(Lf*ls))));
}

double Dynein_bothbound::get_fby() {
  return 0;
}


/*** Get Cartesian Velocities ***/
double Dynein_bothbound::get_d_nmx() {
  return -ls * sin(nba) * get_d_nba();
}

double Dynein_bothbound::get_d_tx() {
  return get_d_nmx() + -lt * sin(nma) * get_d_nma();
}

double Dynein_bothbound::get_d_fmx() {
  return -ls * sin(fba) * get_d_fba();
}

double Dynein_bothbound::get_d_nmy() {
  return ls * cos(nba) * get_d_nba();
}

double Dynein_bothbound::get_d_ty() {
  return get_d_nmy() + lt * cos(nma) * get_d_nma();
}

double Dynein_bothbound::get_d_fmy() {
  return ls * cos(fba) * get_d_fba();
}

/*** Get forces ***/
forces Dynein_bothbound::get_internal() {
  return f;
}

forces Dynein_bothbound::get_brownian() {
  return r;
}

/*** Get energies ***/

double Dynein_bothbound::get_PE() {
  return 0.5*cb*square(nba - eq.nba) + 0.5*cb*square(nma - eq.nma)
    + 0.5*cb*square(fma - eq.fma) + 0.5*cb*square(fba - eq.fba);
}

double Dynein_bothbound::get_KE() {
  return 0;
}

void Dynein_bothbound::log(double t, FILE* data_file) {
  fprintf(data_file, "%.2g\t%.2g\t%.2g\t%.5g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f"
	  "\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",
          get_KE(), get_PE(), get_KE() + get_PE(), t, get_nbx(), get_nby(), get_nmx(), get_nmy(),
          get_tx(), get_ty(), get_fmx(), get_fmy(), get_fbx(), get_fby(), state);
}

void Dynein_bothbound::log_run(float runtime, FILE* data_file) {
  // float run_length = (get_nbx() + get_fbx()) / 2;
  // float ave_step_dist = distance_traveled / steps;
  // float ave_step_time = runtime / steps;

  printf("\n\n***********Run data**********\n");
  // printf("Run length: %f nm\n", run_length);
  // printf("Distance traveled: %f nm\n", distance_traveled);
  // printf("Steps: %d\n", steps);
  // printf("Average step length: %f nm\n", ave_step_dist);
  // printf("Average step time: %g s\n\n\n", ave_step_time);
  // fprintf(data_file, "Run length \tDistance traveled \tSteps \tAve step length \tAve step time\n");
  // fprintf(data_file, "%f\t%f\t%d\t%f\t%g\n", run_length, distance_traveled, steps,
  // 	  ave_step_dist, ave_step_time);
  fclose(data_file);
}
