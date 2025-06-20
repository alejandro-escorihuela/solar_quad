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
  int i, Nf, Nmod, tam = np*3*2, punts = 100;
  quad H0, HF, erH;
  quad e[tam], ec[tam], zb[tam], zc[tam], zp[tam], x_adj[tam], y_adj[tam];
  
  *maxerH = 0.0;
  zeros(e, tam);
  adjunt(x_adj, x, sp);
  adjunt(y_adj, y, sp);

  H0 = hH(z, params, np);
  cart2jacobi(z, zb, params, np);
  pasAB(zb, e, params, np, h, x, y, sp, pA, pB);
  Nf = Nm - 2*abs((int) round(suma(x, sp)));
  Nmod = Nf > punts ? Nf/punts : 1;
  for (i = 0; i < Nf; i++) {
    pasABA(zb, e, params, np, h, a, b, s, pA, pB);
    if (i % Nmod == 0) {
      copy(zc, zb, tam);
      copy(ec, e, tam);
      pasBA(zc, ec, params, np, h, x_adj, y_adj, sp, pA, pB);
      jacobi2cart(zp, zc, params, np);
      HF = hH(zp, params, np);
      erH = fabsq((HF - H0)/H0);
      if (erH > *maxerH)
	*maxerH = erH;
    }
  }
  pasBA(zb, e, params, np, h, x_adj, y_adj, sp, pA, pB);
  jacobi2cart(z, zb, params, np);
}

void evolABAsolar_errH(const char ** z, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerH) {
  int tam = np*3*2;
  quad zq[tam], paramsq[5*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  quad err;

  hq = real2quad(h);
  str2quadV(zq, z, tam);
  str2quadV(paramsq, params, 5*np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);
  evolABAsc(zq, paramsq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &err);
  *lerH = quad2real(log10q(err));

  /* quad2strV(z, zq, tam); */
}

void evolABAsolar_errQ(const char ** z, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerQ) {  
  real meh;
  evolABAsolar_errHQ(z, params, np, Nm, h, a, b, s, x, y, sp, &meh, lerQ);
}

void evolABAsolar_errHQ(const char ** z, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerH, real * lerQ) {
  int tam = np*3*2;
  quad meh, meq, zt[tam];
  quad zq[tam], paramsq[5*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  
  hq = real2quad(h);
  str2quadV(zq, z, tam);
  str2quadV(paramsq, params, 5*np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);
  
  copy(zt, zq, tam);
  evolABAsc(zq, paramsq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh);
  *lerH = quad2real(log10q(meh));
  evolABAsc(zt, paramsq, np, 2*(Nm - abs((int) round(suma(xq, sp)))), hq/2, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh);
  meq = difnorm2(zq, zt, tam/2);
  *lerQ = quad2real(log10q(meq));
  
  /* quad2strV(z, zq, tam); */
}
