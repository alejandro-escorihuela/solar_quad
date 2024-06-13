/* 10-06-2024 */
/* alex */
/* evol.c */
#include <stdio.h>
#include <stdlib.h>
#include "evol.h"

void psc(fluxe phi, quad * z, quad * e, quad * params, quad h, int np) {
  int i, tam = np*3*2;
  quad a[tam], d[tam];

  zeros(d, tam);
  copy(a, z, tam);
  phi(z, d, params, h, np); 
  for (i = 0; i < tam; i++) {
    e[i] += d[i];
    z[i] = a[i] + e[i];
    e[i] += (a[i] - z[i]);
  }

}

void pasABA(quad * z, quad * e, quad * params, int np, quad h, quad * a, quad * b, int s, fluxe pA, fluxe pB) {
  int j;
  for (j = 0; j < s; j++) {
    psc(pA, z, e, params, a[j]*h, np);
    psc(pB, z, e, params, b[j]*h, np);
  }
  psc(pA, z, e, params, a[s]*h, np);
}

void pasAB(quad * z, quad * e, quad * params, int np, quad h, quad * x, quad * y, int sp, fluxe pA, fluxe pB) {
  int j;
  for (j = 0; j < sp; j++) {
    psc(pA, z, e, params, x[j]*h, np);
    psc(pB, z, e, params, y[j]*h, np);
  }
}

void pasBA(quad * z, quad * e, quad * params, int np, quad h, quad * x, quad * y, int sp, fluxe pA, fluxe pB) {
  int j;
  for (j = 0; j < sp; j++) {
    psc(pB, z, e, params, y[j]*h, np);
    psc(pA, z, e, params, x[j]*h, np);
  }
}

void evolABAsc(quad * z, quad * params, int np, int Nm, quad h, quad * a, quad * b, int s, quad * x, quad * y, int sp, fluxe pA, fluxe pB, qcons hH, quad * maxerH) {
  int i, Nf, tam = np*3*2;
  quad H0, HF, erH;
  quad e[tam], ec[tam], zc[tam], x_adj[tam], y_adj[tam];
  
  *maxerH = 0.0;
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
    erH = fabsq((HF - H0)/H0);
    if (erH > *maxerH)
       *maxerH = erH;
  }
  pasBA(z, e, params, np, h, x_adj, y_adj, sp, pA, pB);
}

void evolABAsolar_errH(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * lerH) {
  int tam = np*3*2;
  quad zq[tam], paramsq[np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  quad err;

  hq = real2quad(h);
  real2quadV(zq, z, tam);
  real2quadV(paramsq, params, np);
  real2quadV(aq, a, s + 1);
  real2quadV(bq, b, s);
  real2quadV(xq, x, sp);
  real2quadV(yq, y, sp);
  
  evolABAsc(zq, paramsq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &err);
  *lerH = quad2real(log10q(err));
  
  quad2realV(z, zq, tam);
}

void evolABAsolar_errQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * lerQ) {
  real meh;
  evolABAsolar_errHQ(z, params, np, Nm, h, a, b, s, x, y, sp, &meh, lerQ);
}

void evolABAsolar_errHQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * lerH, real * lerQ) {
  int tam = np*3*2;
  quad meh, meq, zt[tam];
  quad zq[tam], paramsq[np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  
  hq = real2quad(h);
  real2quadV(zq, z, tam);
  real2quadV(paramsq, params, np);
  real2quadV(aq, a, s + 1);
  real2quadV(bq, b, s);
  real2quadV(xq, x, sp);
  real2quadV(yq, y, sp);
  
  copy(zt, zq, tam);
  evolABAsc(zq, paramsq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh);
  *lerH = quad2real(log10q(meh));
  evolABAsc(zt, paramsq, np, 2*(Nm - abs((int) round(suma(xq, sp)))), hq/2, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh);
  meq = difnorm2(zq, zt, tam/2);
  *lerQ = quad2real(log10q(meq));
  
  quad2realV(z, zq, tam);
}
