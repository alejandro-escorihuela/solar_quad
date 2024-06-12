/* 10-06-2024 */
/* alex */
/* evol.c */
#include <stdio.h>
#include <stdlib.h>
#include "evol.h"

void psc(fluxe phi, real * z, real * e, real * params, real h, int np) {
  int i, tam = np*3*2;
  real a[tam], d[tam];

  zeros(d, tam);
  copy(a, z, tam);
  phi(z, d, params, h, np); 
  for (i = 0; i < tam; i++) {
    e[i] += d[i];
    z[i] = a[i] + e[i];
    e[i] += (a[i] - z[i]);
  }

}

void pasABA(real * z, real * e, real * params, int np, real h, real * a, real * b, int s, fluxe pA, fluxe pB) {
  int j;
  for (j = 0; j < s; j++) {
    psc(pA, z, e, params, a[j]*h, np);
    psc(pB, z, e, params, b[j]*h, np);
  }
  psc(pA, z, e, params, a[s]*h, np);
}

void pasAB(real * z, real * e, real * params, int np, real h, real * x, real * y, int sp, fluxe pA, fluxe pB) {
  int j;
  for (j = 0; j < sp; j++) {
    psc(pA, z, e, params, x[j]*h, np);
    psc(pB, z, e, params, y[j]*h, np);
  }
}

void pasBA(real * z, real * e, real * params, int np, real h, real * x, real * y, int sp, fluxe pA, fluxe pB) {
  int j;
  for (j = 0; j < sp; j++) {
    psc(pB, z, e, params, y[j]*h, np);
    psc(pA, z, e, params, x[j]*h, np);
  }
}

void evolABAsc(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, fluxe pA, fluxe pB, qcons hH, real * maxerH) {
  int i, Nf, tam = np*3*2;
  real H0, HF, erH;
  real e[tam], ec[tam], zc[tam], x_adj[tam], y_adj[tam];
  
  *maxerH = 0.0L;
  zeros(e, tam);
  adjunt(x_adj, x, sp);
  adjunt(y_adj, y, sp);

  H0 = hH(z, params, np);
  pasAB(z, e, params, np, h, x, y, sp, pA, pB);
  Nf = Nm - 2*abs((int) round(suma(x, sp)));
  for (i = 0; i < Nf; i++) {
    pasABA(z, e, params, np, h, a, b, s, pA, pB);
    copy(zc, z, tam);
    copy(ec, e, tam);
    pasBA(zc, ec, params, np, h, x_adj, y_adj, sp, pA, pB);
    HF = hH(zc, params, np);
    erH = fabsl((HF - H0)/H0);
    if (erH > *maxerH)
       *maxerH = erH;
  }
  pasBA(z, e, params, np, h, x_adj, y_adj, sp, pA, pB);
}

void evolABAsolar_errH(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * erH) {
  evolABAsc(z, params, np, Nm, h, a, b, s, x, y, sp, phiscHK, phiscHI, ham, erH);
}

void evolABAsolar_errQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * erQ) {
  real meh;
  evolABAsolar_errHQ(z, params, np, Nm, h, a, b, s, x, y, sp, &meh, erQ);
}

void evolABAsolar_errHQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * erH, real * erQ) {
  int tam = np*3*2;
  real meh, zt[tam];
  
  copy(zt, z, tam);
  evolABAsc(z, params, np, Nm, h, a, b, s, x, y, sp, phiscHK, phiscHI, ham, erH);
  evolABAsc(zt, params, np, 2*(Nm - abs((int) round(suma(x, sp)))), h/2, a, b, s, x, y, sp, phiscHK, phiscHI, ham, &meh);
  *erQ = difnorm2(z, zt, tam/2);
}
