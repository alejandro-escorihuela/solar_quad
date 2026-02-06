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

void evolABAsc(quad * z, int nz, quad * params, int np, int Nm, quad h, quad * a, quad * b, int s, quad * x, quad * y, int sp, fluxe pA, fluxe pB, qcons hH, quad * maxerH, const char * nom_arxiu, int p_impr) {
  int i, j, Nf, Nmod, tam = np*3*2, punts;
  quad H0, HF, erH, t, sumcP;
  quad e[tam], ec[tam], zc[tam], zp[tam], x_adj[tam], y_adj[tam];
  FILE *arxiu;

  if (nom_arxiu != NULL) {
    arxiu = fopen(nom_arxiu, "w");
    if (arxiu == NULL) {
      printf("No s'ha pogut crear/sobrescriure l'arxiu %s\n", nom_arxiu);
      exit(-1);
    }
  }
  
  *maxerH = 0.0;
  zeros(e, tam);
  adjunt(x_adj, x, sp);
  adjunt(y_adj, y, sp);

  H0 = hH(z, params, np);
  pasAB(z, e, params, np, h, x, y, sp, pA, pB);
  sumcP = round(suma(x, sp));
  t = h*sumcP;
  Nf = Nm - 2*abs((int) sumcP);
  punts = p_impr > 0 ? p_impr : 100;
  Nmod = Nf > punts ? Nf/punts : 1;
  /* printf("h = %.34Qe,  t = %.34Qe, sumcP = %.34Qe\n", h, t, sumcP); */
  /* exit(-1); */
  for (i = 0; i < Nf; i++) {
    pasABA(z, e, params, np, h, a, b, s, pA, pB);
    t += h;
    if (i % Nmod == 0) {
      copy(zc, z, tam);
      copy(ec, e, tam);
      pasBA(zc, ec, params, np, h, x_adj, y_adj, sp, pA, pB);
      HF = hH(zc, params, np);
      erH = fabsq((HF - H0)/H0);
      if (nom_arxiu != NULL) {
	fprintf(arxiu, "%.34Qe %.34Qe %.34Qe %.34Qe ", t + h*sumcP, H0, HF, erH);
	jacobi2cart(zp, zc, params, np); // aquesta línia hi hauria que generalitzar-la
	for (j = 0; j < tam - 1; j++)
	  fprintf(arxiu, "%.34Qe ", zp[j]);
	fprintf(arxiu, "%.34Qe\n", zp[tam - 1]);
      }
      if (erH > *maxerH)
	*maxerH = erH;
    }
  }
  pasBA(z, e, params, np, h, x_adj, y_adj, sp, pA, pB);

  if (nom_arxiu != NULL)
    fclose(arxiu);
}

void evolABAsolar(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, const char * nom_arxiu, int p_impr) {
  quad zq[nz], zbq[nz], paramsq[np], masq[3*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  quad err;

  hq = real2quad(h);
  str2quadV(zq, z, nz);
  str2quadV(paramsq, params, np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);

  expand_masses(masq, paramsq, np);
  centrar(zq, masq, np);
  cart2jacobi(zbq, zq, masq, np);
  evolABAsc(zbq, nz, masq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &err, nom_arxiu, p_impr);
  jacobi2cart(zq, zbq, masq, np);
}

void evolABAsolar_errH(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerH) {
  quad zq[nz], zbq[nz], paramsq[np], masq[3*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  quad err;

  hq = real2quad(h);
  str2quadV(zq, z, nz);
  str2quadV(paramsq, params, np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);

  expand_masses(masq, paramsq, np);
  centrar(zq, masq, np);
  cart2jacobi(zbq, zq, masq, np);
  evolABAsc(zbq, nz, masq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &err, NULL, 0);
  jacobi2cart(zq, zbq, masq, np);
  *lerH = quad2real(log10q(err));

  /* quad2strV(z, zq, nz); */
}

void evolABAsolar_errQ(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerQ) {  
  real meh;
  evolABAsolar_errHQ(z, nz, params, np, Nm, h, a, b, s, x, y, sp, &meh, lerQ);
}

void evolABAsolar_errHQ(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerH, real * lerQ) {
  quad meh, meq;
  quad zq[nz], zt[nz], zbq[nz], zbt[nz], paramsq[np], masq[3*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  
  hq = real2quad(h);
  str2quadV(zq, z, nz);
  str2quadV(paramsq, params, np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);

  expand_masses(masq, paramsq, np);
  centrar(zq, masq, np);
  cart2jacobi(zbq, zq, masq, np);
  copy(zbt, zbq, nz);
  evolABAsc(zbq, nz, masq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh, NULL, 0);
  *lerH = quad2real(log10q(meh));
  evolABAsc(zbt, nz, masq, np, 2*(Nm - abs((int) round(suma(xq, sp)))), hq/2, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh, NULL, 0);
  jacobi2cart(zq, zbq, masq, np);
  jacobi2cart(zt, zbt, masq, np);
  meq = difnorm2(zq, zt, nz/2);
  *lerQ = quad2real(log10q(meq));
  
  /* quad2strV(z, zq, nz); */
}

void evolABAkpert(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, const char * nom_arxiu, int p_impr) {
  quad zq[nz], zbq[nz], paramsq[np], masq[3*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  quad err;

  hq = real2quad(h);
  str2quadV(zq, z, nz);
  str2quadV(paramsq, params, np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);

  evolABAsc(zbq, nz, masq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscH0_kp, phiscH1_kp, ham, &err, nom_arxiu, p_impr);
}

void evolABAkpert_errH(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerH) {
  quad zq[nz], zbq[nz], paramsq[np], masq[3*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  quad err;

  hq = real2quad(h);
  str2quadV(zq, z, nz);
  str2quadV(paramsq, params, np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);

  evolABAsc(zbq, nz, masq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscH0_kp, phiscH1_kp, ham, &err, NULL, 0);
  *lerH = quad2real(log10q(err));
}

void evolABAkpert_errQ(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerQ) {  
  real meh;
  evolABAkpert_errHQ(z, nz, params, np, Nm, h, a, b, s, x, y, sp, &meh, lerQ);
}

void evolABAkpert_errHQ(const char ** z, int nz, const char ** params, int np, int Nm, real h, const char ** a, const char ** b, int s, const char ** x, const char ** y, int sp, real * lerH, real * lerQ) {
  quad meh, meq;
  quad zq[nz], zt[nz], zbq[nz], zbt[nz], paramsq[np], masq[3*np], hq, aq[s + 1], bq[s], xq[sp], yq[sp];
  
  hq = real2quad(h);
  str2quadV(zq, z, nz);
  str2quadV(paramsq, params, np);
  str2quadV(aq, a, s + 1);
  str2quadV(bq, b, s);
  str2quadV(xq, x, sp);
  str2quadV(yq, y, sp);

  copy(zt, zq, nz);
  evolABAsc(zq, nz, masq, np, Nm, hq, aq, bq, s, xq, yq, sp, phiscH0_kp, phiscH1_kp, ham, &meh, NULL, 0);
  *lerH = quad2real(log10q(meh));
  evolABAsc(zt, nz, masq, np, 2*(Nm - abs((int) round(suma(xq, sp)))), hq/2, aq, bq, s, xq, yq, sp, phiscHK, phiscHI, ham, &meh, NULL, 0);
  meq = difnorm2(zq, zt, nz/2);
  *lerQ = quad2real(log10q(meq));
}
