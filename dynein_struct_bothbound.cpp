#include <stdlib.h>
#include <fstream>
#include <cassert>

#include "dynein_struct_bothbound.h"
#include "dynein_data.h"

/* ********************* BOTHBOUND DYNEIN FUNCTIONS ****************************** */

Dynein::Dynein_bothbound(double nma_init, double fma_init, double nbx_init, double nby_init,
	       bothbound_forces* internal_test, bothbound_forces* brownian_test,
	       bothbound_equilibrium_angles* eq_angles) {
  nbx = nbx_init;
  nby = nby_init;

  nma = nma_init;
  fma = fma_init;

  internal_testcase = internal_test;
  brownian_testcase = brownian_test;
  
  if (eq_angles) {
    eq = *eq_angles; // use test angles
  } else {
    eq = bothbound_pre_powerstroke_internal_angles; // use experimental angles
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

void Dynein_bothbound::set_state(State s) {
  state = s;
}

void Dynein_bothbound::update_velocities() {
  update_internal_forces();
  update_brownian_forces();

  if (state != BOTHBOUND) {
    update_velocities_onebound();
  } else if (state == BOTHBOUND) {
    update_velocities_bothbound();
  }
}

void Dynein_bothbound::switch_to_bothbound() {
  double temp_nma;
  double temp_fma;
  
  if (state == NEARBOUND) {
    temp_nma = M_PI + bma - bba;
    temp_fma = M_PI + uma - uba;
    nbx = get_bbx();
    L = fbx - nbx;
  } else {
    temp_nma = M_PI + uma - uba;
    temp_fma = M_PI + bma - bba;
    nbx = get_fbx();
    L = bbx - fbx;
  }

  distance_traveled += fabs(get_ubx() - get_bbx());
  steps++;
  
  state = BOTHBOUND;
  nma = temp_nma;
  fma = temp_fma;
}

void Dynein_bothbound::switch_to_nearbound() {
  //nearbound -> bma is nma
  double temp_bba = get_nba();
  double temp_bma = nma;
  double temp_uma = fma;
  double temp_uba = get_fba();

  bbx = get_nbx();
  bby = 0;
  
  bba = temp_bba;
  bma = temp_bma;
  uma = temp_uma;
  uba = temp_uba;

  state = NEARBOUND;
}

void Dynein_bothbound::switch_to_farbound() {
  //nearbound -> bma is fma
  double temp_bba = get_fba();
  double temp_bma = fma;
  double temp_uma = nma;
  double temp_uba = get_nba();

  bbx = get_fbx();
  bby = 0;
  
  bba = temp_bba;
  bma = temp_bma;
  uma = temp_uma;
  uba = temp_uba;

  state = FARBOUND;
}

void Dynein_bothbound::unbind() {
  state = UNBOUND;
}

void Dynein_bothbound::update_velocities_onebound() {
  assert(state != BOTHBOUND);
  float A1, A2, A3, A4;  // Start solving for velocities with matrix solution in derivation.pdf
  float B1, B2, B3, B4;
  float C1, C2, C3, C4;
  float D1, D2, D3, D4;
  float X1, X2, X3, X4;
  float Nbb, Nml, Nmr, Nbr;
  float D;
  
  A1 = -4*ls;
  A2 = -3*lt*(sin(bma)*sin(bba) + cos(bma)*cos(bba));
  A3 = 2*lt*(sin(uma)*sin(bba) + cos(uma)*cos(bba));
  A4 = ls*(sin(uba)*sin(bba) + cos(uba)*cos(bba));
  B1 = -3*ls*(sin(bba)*sin(bma) + cos(bba)*cos(bma));
  B2 = -3*lt;
  B3 = +2*lt*(sin(uma)*sin(bma) + cos(uma)*cos(bma));
  B4 = +ls*(sin(uba)*sin(bma) + cos(uba)*cos(bma));
  C1 = -2*ls*(sin(bba)*sin(uma) + cos(bba)*cos(uma));
  C2 = -2*lt*(sin(bma)*sin(uma) + cos(bma)*cos(uma));
  C3 = 2*lt;
  C4 = ls*(sin(uba)*sin(uma) + cos(uba)*cos(uma));
  D1 = ls*(cos(uba)*cos(bba) + sin(uba)*sin(bba));
  D2 = lt*(cos(uba)*cos(bma) + sin(uba)*sin(bma));
  D3 = -lt*(cos(uba)*cos(uma) + sin(uba)*sin(uma));
  D4 = -ls;
  
  X1 = (- 1/gm*f.bmy - 1/gt*f.ty - 1/gm*f.umy - 1/gb*f.uby - r.bmy - r.ty - r.umy - r.uby)*cos(bba)
      + ( 1/gm*f.bmx + 1/gt*f.tx + 1/gm*f.umx + 1/gb*f.ubx + r.bmx + r.tx + r.umx + r.ubx )*sin(bba);
  
  X2 = (-1/gt*f.ty - 1/gm*f.umy - 1/gb*f.uby - r.ty - r.umy - r.uby)*cos(bma) + (1/gt*f.tx + 1/gm*f.umx + 1/gb*f.ubx + r.tx + r.umx + r.ubx)*sin(bma);

  X3 = (-r.umy -r.uby - 1/gm*f.umy - 1/gb*f.uby)*cos(uma) + (r.umx + r.ubx + 1/gm*f.umx + 1/gb*f.ubx)*sin(uma);
  
  X4 = (r.uby + 1/gb*f.uby)*cos(uba) - (r.ubx + 1/gb*f.ubx)*sin(uba);

  Nbb = (-B2*C4*D3*X1 + B2*C3*D4*X1 + A4*C3*D2*X2 - A3*C4*D2*X2 - A4*C2*D3*X2 + A2*C4*D3*X2 + A3*C2*D4*X2 - A2*C3*D4*X2 + A4*B2*D3*X3 - A3*B2*D4*X3 - A4*B2*C3*X4 + A3*B2*C4*X4 + B4*(-C3*D2*X1 + C2*D3*X1 + A3*D2*X3 - A2*D3*X3 - A3*C2*X4 + A2*C3*X4) + B3*(C4*D2*X1 - C2*D4*X1 - A4*D2*X3 + A2*D4*X3 + A4*C2*X4 - A2*C4*X4));

  Nml = (B1*C4*D3*X1 - B1*C3*D4*X1 - A4*C3*D1*X2 + A3*C4*D1*X2 + A4*C1*D3*X2 - A1*C4*D3*X2 - A3*C1*D4*X2 + A1*C3*D4*X2 - A4*B1*D3*X3 + A3*B1*D4*X3 + A4*B1*C3*X4 - A3*B1*C4*X4 + B4*(C3*D1*X1 - C1*D3*X1 - A3*D1*X3 + A1*D3*X3 + A3*C1*X4 - A1*C3*X4) + B3*(-C4*D1*X1 + C1*D4*X1 + A4*D1*X3 - A1*D4*X3 - A4*C1*X4 + A1*C4*X4));

  Nmr = (-B1*C4*D2*X1 + B1*C2*D4*X1 + A4*C2*D1*X2 - A2*C4*D1*X2 - A4*C1*D2*X2 + A1*C4*D2*X2 + A2*C1*D4*X2 - A1*C2*D4*X2 + A4*B1*D2*X3 - A2*B1*D4*X3 - A4*B1*C2*X4 + A2*B1*C4*X4 + B4*(-C2*D1*X1 + C1*D2*X1 + A2*D1*X3 - A1*D2*X3 - A2*C1*X4 + A1*C2*X4) + B2*(C4*D1*X1 - C1*D4*X1 - A4*D1*X3 + A1*D4*X3 + A4*C1*X4 - A1*C4*X4));

  Nbr = (B1*C3*D2*X1 - B1*C2*D3*X1 - A3*C2*D1*X2 + A2*C3*D1*X2 + A3*C1*D2*X2 - A1*C3*D2*X2 - A2*C1*D3*X2 + A1*C2*D3*X2 - A3*B1*D2*X3 + A2*B1*D3*X3 + A3*B1*C2*X4 - A2*B1*C3*X4 + B3*(C2*D1*X1 - C1*D2*X1 - A2*D1*X3 + A1*D2*X3 + A2*C1*X4 - A1*C2*X4) + B2*(-C3*D1*X1 + C1*D3*X1 + A3*D1*X3 - A1*D3*X3 - A3*C1*X4 + A1*C3*X4));

  D = A2*B4*C3*D1 - A2*B3*C4*D1 - A1*B4*C3*D2 + A1*B3*C4*D2 - A2*B4*C1*D3 + A1*B4*C2*D3 + A2*B1*C4*D3 - A1*B2*C4*D3 + A4*(B3*C2*D1 - B2*C3*D1 - B3*C1*D2 + B1*C3*D2 + B2*C1*D3 - B1*C2*D3)+ A2*B3*C1*D4 - A1*B3*C2*D4 - A2*B1*C3*D4 + A1*B2*C3*D4 + A3*(-B4*C2*D1 + B2*C4*D1 + B4*C1*D2 - B1*C4*D2 - B2*C1*D4 + B1*C2*D4);
  
  assert(D != 0);

  d_bba = Nbb/D;
  d_bma = Nml/D;
  d_uma = Nmr/D;
  d_uba = Nbr/D;
}

void Dynein_bothbound::update_velocities_bothbound() {

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
  
  double Cxn = - pm * ls\frac{\sqrt{1 - \frac{ls^4 + lt^4 - 2ls^2lt^2}{4ls^2lf^2} - \frac{lf^2 + 2ls^2 - 2lt^2}{4ls^2}}}{2\sqrt{1 - \frac{L^2}{4lf^2} - \frac{ln^4 - 2L^2ln^2}{4L^2lf^2} - \frac{lf^2 + 2L^2 - 2ln^2}{4L^2}}}\left(-\frac{ln^3 - L^2ln}{L^2lf^2} + \frac{ln}{L^2}\right);
  
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

  d_ln = 1; // from Mathematica
  d_lf = 1; // from Mathematica
}

double Dynein_bothbound::get_binding_rate() {
  if (get_uby() < MICROTUBULE_BINDING_DISTANCE) {
    double bound_energy = 0.5*cb*square(fba - eq.bba);
    return 1e10*exp(-bound_energy/kb/T); // per second
  } else {
    return 0;
  }
}

double Dynein_bothbound::get_unbinding_rate() {
  if (f.bby + r.bby >= UNBINDING_FORCE) { // bad, doesn't take into account forces on other domains?
    printf("unbinding...\n");
    return 1.0;
  } else return 0.0;
}

/*** Set positions, velocities and forces ***/
void Dynein_bothbound::set_bbx(double d) {   // onebound
  assert(C != BOTHBOUND);
  bbx = d;
}

void Dynein_bothbound::set_bby(double d) {   // onebound
  assert(C != BOTHBOUND);
  bby = d;
}

void Dynein_bothbound::set_bba(double d) {   // onebound
  assert(C != BOTHBOUND);
  bba = d;
}

void Dynein_bothbound::set_bma(double d) {   // onebound
  assert(C != BOTHBOUND);
  bma = d;
}

void Dynein_bothbound::set_fma(double d) {   // onebound
  assert(C != BOTHBOUND);
  fma = d;
}

void Dynein_bothbound::set_fba(double d) {   // onebound
  assert(C != BOTHBOUND);
  fba = d;
}

void Dynein_bothbound::set_nma(double d) {   // bothbound
  assert(C == BOTHBOUND);
  nma = d;
}

void Dynein_bothbound::set_fma(double d) {   // bothbound
  assert(C == BOTHBOUND);
  fma = d;
}

void Dynein_bothbound::set_L(double d) {     // bothbound
  assert(C == BOTHBOUND);
  L = d;
}

/*** Angular Velocities ***/

double Dynein_bothbound::get_d_bba() {
  assert(C != BOTHBOUND);
  return d_bba;
}

double Dynein_bothbound::get_d_bma() {
  assert(C != BOTHBOUND);
  return d_bma;
}

double Dynein_bothbound::get_d_fma() {
  assert(C != BOTHBOUND);
  return d_fma;
}

double Dynein_bothbound::get_d_fba() {
  assert(C != BOTHBOUND);
  return d_fba;
}

double Dynein_bothbound::get_d_nma() {  // bothbound
  assert(C == BOTHBOUND);
  int pm = (nma > M_PI) ? -1 : 1; // sign of d_nma depends on value of nma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - ln*ln) / (2*lt*ls),2))
    * (ln / (lt*ls)) * d_ln;
}

double Dynein_bothbound::get_d_fma() {  // bothbound
  assert(C == BOTHBOUND);
  int pm = (fma > M_PI) ? -1 : 1; // sign of d_fma depends on value of fma
  return pm * 1 / sqrt(1 - pow((lt*lt + ls*ls - lf*lf) / (2*lt*ls),2))
    * (lf / (lt*ls)) * d_lf;
}

/*** Get coordinates ***/

double Dynein_bothbound::get_bbx() {
  assert(C != BOTHBOUND);
  return bbx;
}

double Dynein_bothbound::get_bby(){
  assert(C != BOTHBOUND);
  return bby;
}

double Dynein_bothbound::get_bmx() {
  assert(C != BOTHBOUND);
  return ls * cos(get_bba()) + bbx;
}

double Dynein_bothbound::get_bmy(){
  assert(C != BOTHBOUND);
  return ls * sin(get_bba()) + bby;
}

double Dynein_bothbound::get_tx() {
  assert(C != BOTHBOUND);
  return ls * cos(get_bba()) + lt * cos(get_bma()) + bbx;
}

double Dynein_bothbound::get_ty(){
  assert(C != BOTHBOUND);
  return ls * sin(get_bba()) + lt * sin(get_bma()) + bby;
}

double Dynein_bothbound::get_fmx() {
  assert(C != BOTHBOUND);
  return ls * cos(get_bba()) + lt * cos(get_bma()) - lt * cos(get_fma()) + bbx;
}

double Dynein_bothbound::get_umy(){
  assert(C != BOTHBOUND);
  return ls * sin(get_bba()) + lt * sin(get_bma()) - lt * sin(get_fma()) + bby;
}

double Dynein_bothbound::get_fbx() {
  assert(C != BOTHBOUND);
  return ls * cos(get_bba()) + lt * cos(get_bma()) - lt * cos(get_fma()) - ls * cos(get_fba()) + bbx;
}

double Dynein_bothbound::get_uby(){
  assert(C != BOTHBOUND);
  return ls * sin(get_bba()) + lt * sin(get_bma()) - lt * sin(get_fma()) - ls * sin(get_fba()) + bby;
}

double Dynein_bothbound::get_nbx() {   // bothbound
  assert(C == BOTHBOUND);
  return nbx;
}

double Dynein_bothbound::get_nmx() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*(L*ln))) p
         + pm*acos((pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) / (2*(ln*ls))));
}

double Dynein_bothbound::get_tx() {   // bothbound
  assert(C == BOTHBOUND);
  return (pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*L);
}

double Dynein_bothbound::get_fmx() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*cos(acos((pow(L, 2) + pow(lf, 2) - pow(ln, 2)) / (2*(L*lf))) 
         + pm*acos((pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) / (2*(lf*ls))));;
}

double Dynein_bothbound::get_fbx() {   // bothbound
  assert(C == BOTHBOUND);
  return L;
}

double Dynein_bothbound::get_nby() {   // bothbound
  assert(C == BOTHBOUND);
  return 0;
}

double Dynein_bothbound::get_nmy() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (nma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(ln, 2) - pow(lf, 2)) / (2*(L*ln))) 
         + pm*acos((pow(ls, 2) + pow(ln, 2) - pow(lt, 2)) / (2*(ln*ls))));
}

double Dynein_bothbound::get_ty() {   // bothbound
  assert(C == BOTHBOUND);
  return ln*sqrt(1 - pow(pow(L, 2) + pow(ln, 2) - pow(lf, 2), 2) 
  / (4*(pow(L, 2)*pow(ln, 2))));
}

double Dynein_bothbound::get_fmy() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (fma > M_PI) ? -1 : 1;
  return ls*sin(acos((pow(L, 2) + pow(lf, 2) - pow(ln, 2)) / (2*(L*lf))) 
         + pm*acos((pow(ls, 2) + pow(lf, 2) - pow(lt, 2)) / (2*(lf*ls))));
}

double Dynein_bothbound::get_fby() {   // bothbound
  assert(C == BOTHBOUND);
  return 0;
}


/*** Get Cartesian Velocities ***/

double Dynein_bothbound::get_d_bbx() {
  assert(C != BOTHBOUND);
  return 0;
}

double Dynein_bothbound::get_d_bmx() {
  assert(C != BOTHBOUND);
  return ls * d_bba * -sin(bba);
}

double Dynein_bothbound::get_d_tx() {
  assert(C != BOTHBOUND);
  return lt * d_bma * -sin(bma) + get_d_bmx();
}

double Dynein_bothbound::get_d_fmx() {
  assert(C != BOTHBOUND);
  return lt * d_fma * sin(fma) + get_d_tx();
}

double Dynein_bothbound::get_d_fbx() {
  assert(C != BOTHBOUND);
  return ls * d_fba * sin(fba) + get_d_fmx();
}

double Dynein_bothbound::get_d_bby() {
  assert(C != BOTHBOUND);
  return 0;
}

double Dynein_bothbound::get_d_bmy() {
  assert(C != BOTHBOUND);
  return ls * d_bba * cos(bba);
}

double Dynein_bothbound::get_d_ty() {
  assert(C != BOTHBOUND);
  return lt * d_bma * cos(bma) + get_d_bmy();
}

double Dynein_bothbound::get_d_umy() {
  assert(C != BOTHBOUND);
  return lt * d_fma * -cos(fma) + get_d_ty();
}

double Dynein_bothbound::get_d_uby() {
  assert(C != BOTHBOUND);
  return ls * d_fba * -cos(fba) + get_d_umy();
}

double Dynein_bothbound::get_d_nbx() {   // bothbound
  assert(C == BOTHBOUND);
  return ;
}

double Dynein_bothbound::get_d_nmx() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (nma > M_PI) ? -1 : 1;
  return ;
}

double Dynein_bothbound::get_d_tx() {   // bothbound
  assert(C == BOTHBOUND);
  return ;
}

double Dynein_bothbound::get_d_fmx() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (fma > M_PI) ? -1 : 1;
  return ;
}

double Dynein_bothbound::get_d_fbx() {   // bothbound
  assert(C == BOTHBOUND);
  return ;
}

double Dynein_bothbound::get_d_nby() {   // bothbound
  assert(C == BOTHBOUND);
  return ;
}

double Dynein_bothbound::get_d_nmy() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (nma > M_PI) ? -1 : 1;
  return ;
}

double Dynein_bothbound::get_d_ty() {   // bothbound
  assert(C == BOTHBOUND);
  return ;
}

double Dynein_bothbound::get_d_fmy() {   // bothbound
  assert(C == BOTHBOUND);
  int pm = (fma > M_PI) ? -1 : 1;
  return ;
}

double Dynein_bothbound::get_d_fby() {   // bothbound
  assert(C == BOTHBOUND);
  return ;
}

/*** Get forces ***/
forces Dynein_bothbound::get_internal() {
  return f;
}

forces Dynein_bothbound::get_brownian() {
  return r;
}

/*** Get angles ***/

double Dynein_bothbound::get_bba() {
  return bba;
}

double Dynein_bothbound::get_bma() {
  return bma;
}

double Dynein_bothbound::get_fma() {
  return fma;
}

double Dynein_bothbound::get_fba() {
  return fba;
}

State Dynein_bothbound::get_state() {
  return state;
}

/*** Get energies ***/

double Dynein_bothbound::get_PE() {
  return 0;
}

double Dynein_bothbound::get_KE() {
  return 0;
}

void Dynein_bothbound::log(double t, FILE* data_file) {
  fprintf(data_file, "%.2g\t%.2g\t%.2g\t%.5g\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\n",
          get_KE(), get_PE(), get_KE() + get_PE(), t, get_bbx(), get_bby(), get_bmx(), get_bmy(),
          get_tx(), get_ty(), get_fmx(), get_umy(), get_fbx(), get_uby(), state);
}

void Dynein_bothbound::log_run(float runtime) {
  FILE* data_file = fopen("run_data.txt", "w");

  float run_length = get_bbx();
  float ave_step_dist = distance_traveled / steps;
  float ave_step_time = runtime / steps;

  printf("\n\n***********Run data**********\n");
  printf("Run length: %f nm\n", run_length);
  printf("Distance traveled: %f nm\n", distance_traveled);
  printf("Steps: %d\n", steps);
  printf("Average step length: %f nm\n", ave_step_dist);
  printf("Average step time: %g s\n\n\n", ave_step_time);
  fprintf(data_file, "Run length \tDistance traveled \tSteps \tAve step length \tAve step time\n");
  fprintf(data_file, "%f\t%f\t%d\t%f\t%g\n", run_length, distance_traveled, steps, ave_step_dist, ave_step_time);
  fclose(data_file);
}

void Dynein_bothbound::resetLog() {
  FILE* data_file = fopen("data.txt", "w");
  FILE* config_file = fopen("config.txt", "w");

  fprintf(config_file, "#gb\tgm\tgt\tdt\truntime?\tstate\n");
  fprintf(config_file, "%g\t%g\t%g\t%g\t%g\t%d\n",
          (double) gb, (double) gm, (double) gt, dt, runtime, (int) state);
  fprintf(data_file,
	  "#KE\tPE\tEnergy\tt\tnbX\tnbY\tnmx\tnmy\ttX\ttY\tfmx\tfmy\tfbx\tuby\tS\n");

  fclose(data_file);
  fclose(config_file);
}
