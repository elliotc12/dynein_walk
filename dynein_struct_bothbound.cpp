#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct_bothbound.h"
#include "dynein_data.h"

/* ********************* BOTHBOUND DYNEIN FUNCTIONS ****************************** */

Dynein::Dynein_bothbound(double nma_init, double fma_init, double nbx_init,
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
  bothbound_forces old_r = old_dynein->get_internal_forces();
  bothbound_forces old_r = old_dynein->get_internal_forces();

  rand = mtrand;

  distance_traveled = old_dynein->distance_traveled;
  steps = old_dynein->steps;

  if (old_dynein->get_state() == State::NEARBOUND) {
    nbx = old_dynein->get_bbx();
    nby = 0;

    f.nbx = old_f.bbx;     f.nby = old_f.bby;
    f.nmx = old_f.bmx;     f.nmy = old_f.bmy;
    f.tx  = old_f.tx;      f.ty  = old_f.ty;
    f.fmx = old_f.umx;     f.fmy = old_f.umy;
    f.fbx = old_f.ubx;     f.fby = old_f.uby;

    r.nbx = old_r.bbx;     r.nby = old_r.bby;
    r.nmx = old_r.bmx;     r.nmy = old_r.bmy;
    r.tx  = old_r.tx;      r.ty  = old_r.ty;
    r.fmx = old_r.umx;     r.fmy = old_r.umy;
    r.fbx = old_r.ubx;     r.fby = old_r.uby;

    nma = M_PI + old_dynein->get_bma() - old_dynein->get_bba();
    fma = M_PI + old_dynein->get_uma() - old_dynein->get_uba();

    L = old_dynein->ubx - old_dynein->bbx;

  } else {
    nbx = old_dynein->get_ubx();
    nby = 0;

    f.nbx = old_f.ubx;     f.nby = old_f.uby;
    f.nmx = old_f.umx;     f.nmy = old_f.umy;
    f.tx  = old_f.tx;      f.ty  = old_f.ty;
    f.fmx = old_f.bmx;     f.fmy = old_f.bmy;
    f.fbx = old_f.bbx;     f.fby = old_f.bby;

    r.nbx = old_r.ubx;     r.nby = old_r.uby;
    r.nmx = old_r.umx;     r.nmy = old_r.umy;
    r.tx  = old_r.tx;      r.ty  = old_r.ty;
    r.fmx = old_r.bmx;     r.fmy = old_r.bmy;
    r.fbx = old_r.bbx;     r.fby = old_r.bby;

    nma = M_PI + old_dynein->get_uma() - old_dynein->get_uba();
    fma = M_PI + old_dynein->get_bma() - old_dynein->get_bba();

    L = old_dynein->bbx - old_dynein->ubx;
  }

  update_velocity();
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

State Dynein_bothbound::get_state() {
  return State::BOTHBOUND;
}

void Dynein_bothbound::update_velocities() {
  update_internal_forces();
  update_brownian_forces();

  int pm_n = 1;
  int pm_f = 1;
  if (nma > M_PI) pm_n = -1;
  if (fma > M_PI) pm_f = -1;

  double L = fbx - nbx;

  double Axn = (pow(L, 2)*pow(lt, 2) - pow(L, 2)*pow(ls, 2)) / (2*(L*pow(ln, 3))) 
  + (lf*pow(lt, 2) - lf*pow(ls, 2)) / (2*(L*pow(ln, 2))) 
  + (pow(lf, 2)*pow(ls, 2) - pow(lf, 2)*pow(lt, 2)) / (2*(L*pow(ln, 3))) 
  - ls*(pm_n*((pow(L, 2) / (2*pow(ln, 3)) 
  + (pow(lf, 4) - 2*(pow(L, 2)*pow(lf, 2))) / (2*(pow(L, 2)*pow(ln, 3))) 
  - ln / (2*pow(L, 2)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(ln, 2)*pow(ls, 2))) 
         + (2*pow(lt, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))) 
      / (2*sqrt(1 - pow(L, 2) / (4*pow(ln, 2)) 
  + (2*(pow(L, 2)*pow(lf, 2)) - pow(lf, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
  + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2)))) 
  - ls*(pm_n*(((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2))) 
  / (2*(pow(ln, 3)*pow(ls, 2))) - ln / (2*pow(ls, 2)))
 *sqrt(1 - pow(L, 2) / (4*pow(ln, 2)) 
         + (2*(pow(L, 2)*pow(lf, 2)) - pow(lf, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
         + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
      / (2*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(ln, 2)*pow(ls, 2))) 
  + (2*pow(lt, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1));

  double Axf = ls*(pm*(((lf*pow(L, 2) - pow(lf, 3)) / (pow(L, 2)*pow(ln, 2)) + lf / pow(L, 2))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(ln, 2)*pow(ls, 2)))
         + (2*pow(lt, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))) 
  / (-2*sqrt(1 - pow(L, 2) / (4*pow(ln, 2)) 
  + (2*(pow(L, 2)*pow(lf, 2)) - pow(lf, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
  + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))));
  
  double Bxn = ln/L;
  
  double Bxf = -lf/L;
  
  double Cxn = - pm * ls\frac{\sqrt{1 - \frac{ls^4 + lt^4 - 2ls^2lt^2}{4ls^2lf^2} -
	\frac{lf^2 + 2ls^2 - 2lt^2}{4ls^2}}}{2\sqrt{1 - \frac{L^2}{4lf^2}
      - \frac{ln^4 - 2L^2ln^2}{4L^2lf^2} - \frac{lf^2 + 2L^2 - 2ln^2}{4L^2}}}
        \left(-\frac{ln^3 - L^2ln}{L^2lf^2} + \frac{ln}{L^2}\right);
  
  double Cxf = (pow(L, 2)*pow(lt, 2) - pow(L, 2)*pow(ls, 2)) / (2*(L*pow(lf, 3))) 
  + (ln*pow(lt, 2) - ln*pow(ls, 2)) / (2*(L*pow(lf, 2))) 
  + (pow(ln, 2)*pow(ls, 2) - pow(ln, 2)*pow(lt, 2)) / (2*(L*pow(lf, 3))) 
  - ls*(pm*((pow(L, 2) / (2*pow(lf, 3)) 
  + (pow(ln, 4) - 2*(pow(L, 2)*pow(ln, 2))) / (2*(pow(L, 2)*pow(lf, 3))) 
  - lf / (2*pow(L, 2)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(lf, 2)*pow(ls, 2))) 
         + (2*pow(lt, 2) - pow(lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))) 
      / (2*sqrt(1 - pow(L, 2) / (4*pow(lf, 2)) 
  + (2*(pow(L, 2)*pow(ln, 2)) - pow(ln, 4)) / (4*(pow(L, 2)*pow(lf, 2))) 
  + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2)))) 
  - ls*(pm*(((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2))) 
  / (2*(pow(lf, 3)*pow(ls, 2))) - lf / (2*pow(ls, 2)))
 *sqrt(1 - pow(L, 2) / (4*pow(lf, 2)) 
         + (2*(pow(L, 2)*pow(ln, 2)) - pow(ln, 4)) / (4*(pow(L, 2)*pow(lf, 2))) 
         + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
      / (2*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(lf, 2)*pow(ls, 2))) 
  + (2*pow(lt, 2) - pow(lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1));
  
  double Ayn = ((pow(lt, 2) - pow(ls, 2) - pow(ln, 2)) / (2*pow(ln, 2)) + 1)
 *sqrt(1 - pow(L, 2) / (4*pow(ln, 2)) 
         + (2*(pow(L, 2)*pow(lf, 2)) - pow(lf, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
         + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))) 
  + (pow(L, 2) / (2*pow(ln, 3)) 
       + (pow(lf, 4) - 2*(pow(L, 2)*pow(lf, 2))) / (2*(pow(L, 2)*pow(ln, 3))) 
       - ln / (2*pow(L, 2)))*(pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) 
      / (4*(ln*sqrt(1 - pow(L, 2) / (4*pow(ln, 2)) 
  + (2*(pow(L, 2)*pow(lf, 2)) - pow(lf, 4)) / (4*pow(L, 2[n] + 2)) 
  + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm*(ls(1 / L + (pow(lf, 2) - pow(L, 2) - pow(ln, 2)) / (2*pow(Lln, 2)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(ln, 2)*pow(ls, 2))) 
         + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + ls*(pm*((pow(L, 2) + pow(ln, 2) - pow(lf, 2))
 *((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2))) 
     / (2*(pow(ln, 3)*pow(ls, 2))) - ln / (2*pow(ls, 2))))) 
      / (4*(Lln
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(ln, 2)*pow(ls, 2))) 
         + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)));
  
  double Ayf = ((lf*pow(L, 2) - pow(lf, 3)) / (pow(L, 2)*pow(ln, 2)) + lf / pow(L, 2))
 *(pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) 
  / (4*(ln*sqrt(1 - pow(L, 2) / (4*pow(ln, 2)) 
  + (2*(pow(L, 2)*pow(lf, 2)) - pow(lf, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
  + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm*(ls(-(lf / (L*ln)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(ln, 2)*pow(ls, 2))) 
         + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + lf*(pm*(pow(L, 2) + pow(ln, 2) - pow(lf, 2))) 
      / (4*(L*(ln*(ls*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(ln, 2)*pow(ls, 2))) 
  + (2*pow(lf, 2) - pow(ln, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)))));
  
  double Byn = (2*ln - ln*(pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / pow(L, 2)) 
  / (2*sqrt(pow(ln, 2) 
  - pow(pow(L, 2) + pow(ln, 2) - pow(lf, 2), 2) / (4*pow(L, 2))));
  
  double Byf = lf*(pow(L, 2) + pow(ln, 2) - pow(lf, 2)) 
  / (2*(pow(L, 2)
 *sqrt(pow(ln, 2) - pow(pow(L, 2) + pow(ln, 2) - pow(lf, 2), 2) / (4*pow(L, 2)))));
  
  double Cyn = ((ln*pow(L, 2) - pow(ln, 3)) / (pow(L, 2)*pow(ln, 2)) + ln / pow(L, 2))
 *(pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) 
  / (4*(lf*sqrt(1 - pow(L, 2) / (4*pow(lf, 2)) 
  + (2*(pow(L, 2)*pow(ln, 2)) - pow(ln, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
  + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm*(ls(-(ln / (L*lf)))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(lf, 2)*pow(ls, 2))) 
         + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + ln*(pm*(pow(L, 2) + pow(lf, 2) - pow(ln, 2))) 
      / (4*(L*(lf*(ls*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(lf, 2)*pow(ls, 2))) 
  + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)))));
  
  double Cyf = ((pow(lt, 2) - pow(ls, 2) - pow(lf, 2)) / (2*pow(lf, 2)) + 1)
 *sqrt(1 - pow(L, 2) / (4*pow(lf, 2)) 
         + (2*(pow(L, 2)*pow(ln, 2)) - pow(ln, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
         + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))) 
  + (pow(L, 2) / (2*pow(lf, 3)) 
       + (pow(ln, 4) - 2*(pow(L, 2)*pow(ln, 2))) / (2*(pow(L, 2)*pow(lf, 3))) 
       - lf / (2*pow(L, 2)))*(pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) 
      / (4*(lf*sqrt(1 - pow(L, 2) / (4*pow(lf, 2)) 
  + (2*(pow(L, 2)*pow(ln, 2)) - pow(ln, 4)) / (4*(pow(L, 2)*pow(ln, 2))) 
  + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(L, 2)) / (4*pow(L, 2))))) 
  + pm*(ls(1 / L + (pow(ln, 2) - pow(L, 2) - pow(lf, 2)) / (2*(L*pow(lf, 2))))
 *sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
         / (4*(pow(lf, 2)*pow(ls, 2))) 
         + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1)) 
  + ls*(pm*((pow(L, 2) + pow(lf, 2) - pow(ln, 2))
 *((pow(ls, 4) + pow(lt, 4) - 2*(pow(ls, 2)*pow(lt, 2)))
     / (2*(pow(lf, 3)*pow(ls, 2))) - lf / (2*pow(ls, 2))))) 
      / (4*(L*(lf*sqrt((2*(pow(ls, 2)*pow(lt, 2)) - pow(ls, 4) - pow(lt, 4)) 
  / (4*(pow(lf, 2)*pow(ls, 2))) 
  + (2*pow(ln, 2) - pow(lf, 2) - 2*pow(ls, 2)) / (4*pow(ls, 2)) + 1))));

  d_ln = (Cyf*Fxmr - Cxf*Fymr + Cyf*g*Rxmr - Cxf*g*Rymr)/((Cxn*Cyf - Cxf*Cyn)*g) - (((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr)) +
    (Cxn*Cyf - Cxf*Cyn)*(-((Byn*Fxmr)/g) + (Cxn*Fyt)/g - Byn*Rxmr + Cxn*Ryt))*(Xbm - Xt))/g -
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn) *(-((Bxn*Fxmr)/g) + (Cxn*Fxt)/g - Bxn*Rxmr +
    Cxn*Rxt) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr))*(Ybm -
    Yt))/g)*(Cyf*Xfm - Cyf*Xt - Cxf*Yfm + Cxf*Yt))/ ((Cxn*Cyf - Cxf*Cyn)*g*(-((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(Ybm - Yt)*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Xfm - Xt))/g) + (Bxn*(-Xfm + Xt))/g) -
    (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) + (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(Xbm - Xt)*((Cxn*Cyf - Cxf*Cyn)* ((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) -
    (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn* (-Yfm + Yt))/g)))/g)) - (((-(Cyf*Xfb) +
    Cyf*Xfm + Cxf*Yfb - Cxf*Yfm)/((Cxn*Cyf - Cxf*Cyn)*g) - (((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm -
    Xt)*((Byn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(Xfb - Xfm))/g -
    (Cxn*(Yfb - Yfm))/g)))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Bxn*(Cxn* Cyf - Cxf*Cyn)*(Xfb - Xfm))/g -
    (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g))*(Ybm - Yt))/g)*(Cyf*Xfm -
    Cyf*Xt - Cxf*Yfm + Cxf*Yt))/((Cxn*Cyf - Cxf*Cyn)*g*(-((Cxn* (Cxn*Cyf - Cxf*Cyn)*(Ybm -
    Yt)*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Xfm - Xt))/g) + (Bxn*(-Xfm + Xt))/ g) - (-(Bxn*Cxf) +
    Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm -
    Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g)))* ((-((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Ayn*Fxmr)/g) + (Cxn*Fyml)/g -
    Ayn*Rxmr + Cxn*Ryml) - (-(Ayn*Cxf) + Ayf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr +
    Cxn*Rymr))*(Xbb - Xbm))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*((Cxn*Fxml)/g -
    (Axn*Fxmr)/g + Cxn*Rxml - Axn*Rxmr) - (-(Axn*Cxf) + Bxf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g -
    Cyn*Rxmr + Cxn*Rymr))* (Ybb - Ybm))/g)*(Ybm - Yt))/g) - (-((-(Byn*Cxf) +
    Byf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr)) + (Cxn*Cyf -
    Cxf*Cyn)*(-((Byn*Fxmr)/g) + (Cxn*Fyt)/g - Byn*Rxmr + Cxn*Ryt))* (-(((Cxn*Cxn)*((Cxn*Cyf -
    Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm + Xt)*(Ybb - Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf -
    Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb - Xbm)*(-Ybm + Yt))/(g*g)))* (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm
    - Yt)*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Xfm - Xt))/g) + (Bxn*(-Xfm + Xt))/g) - (-(Bxn*Cxf) +
    Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm -
    Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) - ((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr)) +
    (Cxn*Cyf - Cxf*Cyn)*(-((Byn*Fxmr)/g) + (Cxn*Fyt)/g - Byn*Rxmr + Cxn*Ryt))* (Xbm - Xt))/g -
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Bxn*Fxmr)/g) + (Cxn*Fxt)/g - Bxn*Rxmr +
    Cxn*Rxt) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr))*(Ybm -
    Yt))/g)* (-((-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm + Xt)*(Ybb -
    Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb - Xbm)*(-Ybm +
    Yt))/(g*g))* ((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g))) - (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm -
    Yt)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbb - Xbm)*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm + Xt))/g -
    (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/ g) + (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(Ybb - Ybm)*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm + Xt))/g - (-(Axn*Cxf) +
    Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g))/g)))/ ((-((Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbb - Xbm)*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g -
    (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g)))/g) + (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*(Ybb - Ybm)*((Axn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Axn*Cxf) +
    Bxf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g)))/g)*(Ybm - Yt))/g) - ((Byn*(Cxn*Cyf -
    Cxf*Cyn)*(Xfb - Xfm))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g))*
    (-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm + Xt)*(Ybb - Ybm))/(g*g)) +
    ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb - Xbm)*(-Ybm + Yt))/(g*g)))*
    (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm - Yt)*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Xfm - Xt))/g) + (Bxn*(-Xfm
    + Xt))/g) - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) +
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm - Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm -
    Yt))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) -
    ((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm - Xt)*((Byn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Byn*Cxf) +
    Byf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g)))/g - (Cxn*(Cxn*Cyf -
    Cxf*Cyn)*((Bxn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(Xfb - Xfm))/g
    - (Cxn*(Yfb - Yfm))/g))*(Ybm - Yt))/g)* (-((-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf -
    Cxf*Cyn))*(-Xbm + Xt)*(Ybb - Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf -
    Cxf*Cyn))*(Xbb - Xbm)*(-Ybm + Yt))/(g*g))* ((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm
    - Yt))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g))) - (Cxn*(Cxn*Cyf
    - Cxf*Cyn)*(Ybm - Yt)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbb - Xbm)*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm +
    Xt))/g - (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/ g) +
    (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybb - Ybm)*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm + Xt))/g - (-(Axn*Cxf) +
    Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g))/g)); // from Mathematica

  
  d_lf = (-(Cyn*Fxmr) + Cxn*Fymr - Cyn*g*Rxmr + Cxn*g*Rymr)/((Cxn*Cyf - Cxf*Cyn)*g) -
           (((Cxn*(Cxn*Cyf - Cxf*Cyn)* (-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g -
           Cyn*Rxmr + Cxn*Rymr)) + (Cxn*Cyf - Cxf*Cyn)*(-((Byn*Fxmr)/g) + (Cxn*Fyt)/g - Byn*Rxmr +
           Cxn*Ryt))*(Xbm - Xt))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Bxn*Fxmr)/g)
           + (Cxn*Fxt)/g - Bxn*Rxmr + Cxn*Rxt) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*Fxmr)/g) +
           (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr))* (Ybm - Yt))/g)*(-(Cyn*Xfm) + Cyn*Xt + Cxn*Yfm -
           Cxn*Yt))/ ((Cxn*Cyf - Cxf*Cyn)*g*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm - Yt)*((Cxn*Cyf -
           Cxf*Cyn)*(-((Cxn*(Xfm - Xt))/g) + (Bxn*(-Xfm + Xt))/g) - (-(Bxn*Cxf) +
           Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm
           - Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) +
           Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g)) - (((Cyn*Xfb - Cyn*Xfm -
           Cxn*Yfb + Cxn*Yfm)/((Cxn*Cyf - Cxf*Cyn)*g) - (((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm - Xt)*
           ((Byn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(Xfb - Xfm))/g -
           (Cxn*(Yfb - Yfm))/g)))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Bxn*(Cxn*Cyf - Cxf*Cyn)*(Xfb -
           Xfm))/g - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g))*(Ybm -
           Yt))/g)* (-(Cyn*Xfm) + Cyn*Xt + Cxn*Yfm - Cxn*Yt))/ ((Cxn*Cyf -
           Cxf*Cyn)*g*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm - Yt)*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Xfm -
           Xt))/g) + (Bxn*(-Xfm + Xt))/g) - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm
           + Yt))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm - Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm +
           Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm +
           Yt))/g)))/g)))* ((-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf -
           Cxf*Cyn)*(-((Ayn*Fxmr)/g) + (Cxn*Fyml)/g - Ayn*Rxmr + Cxn*Ryml) - (-(Ayn*Cxf) +
           Ayf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr))*(Xbb - Xbm))/g) +
           (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*((Cxn*Fxml)/g - (Axn*Fxmr)/g + Cxn*Rxml -
           Axn*Rxmr) - (-(Axn*Cxf) + Bxf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr +
           Cxn*Rymr))* (Ybb - Ybm))/g)*(Ybm - Yt))/g) - (-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*Fxmr)/g) +
           (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr)) + (Cxn*Cyf - Cxf*Cyn)*(-((Byn*Fxmr)/g) + (Cxn*Fyt)/g
           - Byn*Rxmr + Cxn*Ryt))* (-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm +
           Xt)*(Ybb - Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb -
           Xbm)*(-Ybm + Yt))/(g*g)))* (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm - Yt)*((Cxn*Cyf -
           Cxf*Cyn)*(-((Cxn*(Xfm - Xt))/g) + (Bxn*(-Xfm + Xt))/g) - (-(Bxn*Cxf) +
           Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm
           - Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) +
           Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g) - ((Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(-((-(Byn*Cxf) + Byf*Cxn)*(-((Cyn*Fxmr)/g) + (Cxn*Fymr)/g - Cyn*Rxmr +
           Cxn*Rymr)) + (Cxn*Cyf - Cxf*Cyn)*(-((Byn*Fxmr)/g) + (Cxn*Fyt)/g - Byn*Rxmr + Cxn*Ryt))*
           (Xbm - Xt))/g - (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Cxn*Cyf - Cxf*Cyn)*(-((Bxn*Fxmr)/g) +
           (Cxn*Fxt)/g - Bxn*Rxmr + Cxn*Rxt) - (-(Bxn*Cxf) + Bxf*Cxn)*(-((Cyn*Fxmr)/g) +
           (Cxn*Fymr)/g - Cyn*Rxmr + Cxn*Rymr))*(Ybm - Yt))/g)* (-((-(((Cxn*Cxn)*((Cxn*Cyf -
           Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm + Xt)*(Ybb - Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf -
           Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb - Xbm)*(-Ybm + Yt))/(g*g))* ((Cxn*Cyf -
           Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm
           + Xt))/g - (Cxn*(-Yfm + Yt))/g))) - (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm - Yt)*(-((Cxn*(Cxn*Cyf
           - Cxf*Cyn)*(Xbb - Xbm)*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm + Xt))/g - (-(Ayn*Cxf) +
           Ayf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/ g) + (Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(Ybb - Ybm)*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm + Xt))/g - (-(Axn*Cxf) +
           Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g))/g)))/ ((-((Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbb - Xbm)*((Ayn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g
           - (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g)))/g) + (Cxn*(Cxn*Cyf
           - Cxf*Cyn)*(Ybb - Ybm)*((Axn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Axn*Cxf) +
           Bxf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g)))/g)*(Ybm - Yt))/g) - ((Byn*(Cxn*Cyf
           - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb -
           Yfm))/g))* (-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm + Xt)*(Ybb -
           Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb - Xbm)*(-Ybm +
           Yt))/(g*g)))* (-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybm - Yt)*((Cxn*Cyf - Cxf*Cyn)*(-((Cxn*(Xfm -
           Xt))/g) + (Bxn*(-Xfm + Xt))/g) - (-(Bxn*Cxf) + Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm
           + Yt))/g)))/g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm - Xt)*((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm +
           Xt))/g - (Cxn*(Yfm - Yt))/g) - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm +
           Yt))/g)))/g) - ((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbm - Xt)*((Byn*(Cxn*Cyf - Cxf*Cyn)*(Xfb -
           Xfm))/g - (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g)))/g -
           (Cxn*(Cxn*Cyf - Cxf*Cyn)*((Bxn*(Cxn*Cyf - Cxf*Cyn)*(Xfb - Xfm))/g - (-(Bxn*Cxf) +
           Bxf*Cxn)*((Cyn*(Xfb - Xfm))/g - (Cxn*(Yfb - Yfm))/g))*(Ybm - Yt))/g)*
           (-((-(((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(-Xbm + Xt)*(Ybb -
           Ybm))/(g*g)) + ((Cxn*Cxn)*((Cxn*Cyf - Cxf*Cyn)*(Cxn*Cyf - Cxf*Cyn))*(Xbb - Xbm)*(-Ybm +
           Yt))/(g*g))* ((Cxn*Cyf - Cxf*Cyn)*((Byn*(-Xfm + Xt))/g - (Cxn*(Yfm - Yt))/g) -
           (-(Byn*Cxf) + Byf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g))) - (Cxn*(Cxn*Cyf -
           Cxf*Cyn)*(Ybm - Yt)*(-((Cxn*(Cxn*Cyf - Cxf*Cyn)*(Xbb - Xbm)*((Ayn*(Cxn*Cyf -
           Cxf*Cyn)*(-Xfm + Xt))/g - (-(Ayn*Cxf) + Ayf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm +
           Yt))/g)))/ g) + (Cxn*(Cxn*Cyf - Cxf*Cyn)*(Ybb - Ybm)*((Axn*(Cxn*Cyf - Cxf*Cyn)*(-Xfm +
           Xt))/g - (-(Axn*Cxf) + Bxf*Cxn)*((Cyn*(-Xfm + Xt))/g - (Cxn*(-Yfm + Yt))/g)))/g))/g));
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
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - ln*ln) / (2*lt*ls),2))
    * (ln / (lt*ls)) * d_ln;
}

double Dynein_bothbound::get_d_fma() {
  int pm = (fma > M_PI) ? -1 : 1; // sign of d_fma depends on value of fma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - lf*lf) / (2*lt*ls),2))
    * (lf / (lt*ls)) * d_lf;
}

/*** Get coordinates ***/
double Dynein_bothbound::get_nbx() {
  return nbx;
}

double Dynein_bothbound::get_nmx() {
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*(L*ln)))
         + pm*acos((pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) / (2*(ln*ls))));
}

double Dynein_bothbound::get_tx() {
  return (pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*L);
}

double Dynein_bothbound::get_fmx() {
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(lf, 2) - pow(ln, 2)) / (2*(L*lf)))
		+ pm*acos((pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) / (2*(lf*ls))));;
}

double Dynein_bothbound::get_fbx() {
  return L;
}

double Dynein_bothbound::get_nby() {
  return 0;
}

double Dynein_bothbound::get_nmy() {
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*(L*ln)))
		+ pm*acos((pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) / (2*(ln*ls))));
}

double Dynein_bothbound::get_ty() {
  return ln*sqrt(1 - pow(pow(L, 2) + pow(ln, 2) - pow(lf, 2), 2)
		 / (4*(pow(L, 2)*pow(ln, 2))));
}

double Dynein_bothbound::get_fmy() {
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(lf, 2) - pow(ln, 2)) / (2*(L*lf)))
		+ pm*acos((pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) / (2*(lf*ls))));
}

double Dynein_bothbound::get_fby() {
  return 0;
}


/*** Get Cartesian Velocities ***/
double Dynein_bothbound::get_d_nmx() {
  return -ls * sin(nba) * d_nba;
}

double Dynein_bothbound::get_d_tx() {
  return get_d_nmx() + -lt * sin(nma) * d_nma;
}

double Dynein_bothbound::get_d_fmx() {
  return -ls * sin(fba) * d_fba;
}

double Dynein_bothbound::get_d_nmy() {
  return ls * cos(nba) * d_nba;
}

double Dynein_bothbound::get_d_ty() {
  return get_d_nmy() + lt * cos(nma) * d_nma;
}

double Dynein_bothbound::get_d_fmy() {
  return ls * cos(fba) * d_fba;
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
  float run_length = (get_nbx() + get_fbx()) / 2;
  float ave_step_dist = distance_traveled / steps;
  float ave_step_time = runtime / steps;

  printf("\n\n***********Run data**********\n");
  printf("Run length: %f nm\n", run_length);
  printf("Distance traveled: %f nm\n", distance_traveled);
  printf("Steps: %d\n", steps);
  printf("Average step length: %f nm\n", ave_step_dist);
  printf("Average step time: %g s\n\n\n", ave_step_time);
  fprintf(data_file, "Run length \tDistance traveled \tSteps \tAve step length \tAve step time\n");
  fprintf(data_file, "%f\t%f\t%d\t%f\t%g\n", run_length, distance_traveled, steps,
	  ave_step_dist, ave_step_time);
  fclose(data_file);
}
