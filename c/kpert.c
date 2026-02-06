/* 03-02-2026 */
/* alex */
/* kpert.c */
#include <stdio.h>
#include <stdlib.h>
#include "kpert.h"

void kepler_sc_kp(quad * z, quad * dz, quad h) {
  /* Sergio Blanes and Fernando Casas: A Concise Introduction to Geometric Numerical Integrator p[28,29] */
  /* Adapted to compensated summation by Ander Murua */
  quad q[2], v[2];
  quad t, mu = 1.0Q, r0, v02, u, a;
  quad c, s, s2, sig, psi, w, x, dx, adx, adx_ant;
  quad fft, gg, fp, gpt, aux;
  int j;

  for (j = 0; j < 2; j++) {
    q[j] = z[j];
    v[j] = z[j + 2];  
  }
  t = h;
  r0 = sqrtq(dot(q, q, 2));
  v02 = dot(v, v, 2);
  u = dot(q, v, 2);
  a = -mu/(v02 - ((2.0*mu)/r0));
  w = sqrtq(mu/(a*a*a));
  sig = 1.0 - r0/a;
  psi = u/(w*a*a);
  adx = INFINITY;
  x = w*t*a/r0;
  do {
    adx_ant = adx;
    c = cosq(x);
    s = sinq(x);
    s2 = sinq(x/2.0);
    s2 = 2.0*s2*s2;
    dx = (x - sig*s + psi*s2 - w*t)/(1.0 - sig*c + psi*s);
    x = x - dx;
    adx = fabsq(dx);
  } while (adx != 0.0 && adx < adx_ant);
  aux = 1.0 - sig*c + psi*s;
  fft = -s2*a/r0;
  gg = t + ((s - x)/w);
  fp = (-a*w*s)/(aux*r0);
  gpt = -s2/aux;
  for (j = 0; j < 2; j++) {
    dz[j] = fft*q[j] + gg*v[j];
    dz[j + 2] = fp*q[j] + gpt*v[j];
  }
}

void phiscH0_kp(quad * z, quad * dz, int nz, quad * par, quad h, int np) {
  (void) nz;
  (void) par;
  (void) np;
  kepler_sc_kp(z, dz, h);
}

void phiscH1_kp(quad * z, quad * dz, int nz, quad * par, quad h, int np) {
  int j;
  quad q02, q12, r, r7, frc, eps = par[0], alp = par[1];

  q02 = z[0]*z[0];
  q12 = z[1]*z[1];
  r = sqrtq(q02 + q12);
  r7 = powq(r, 7);
  frc = -h*1.5*eps/r7;
  dz[0] = dz[1] = 0.0;
  dz[2] = z[0]*((1 - 3*alp)*q02 + (1 + 2*alp)*q12)*frc;
  dz[3] = z[1]*((1 - 5*alp)*q02 + q12)*frc;
}

quad ham_kp(quad * z, int nz, quad * par, int np) {
  quad q[2], v[2];
  quad eps = par[0], alp = par[1];
  quad r, r2, r3, A, eB;
  int j;
  
  for (j = 0; j < 2; j++) {
    q[j] = z[j];
    v[j] = z[j + 2];
  }
  r = sqrtq(q[0]*q[0] + q[1]*q[1]);
  r2 = r*r;
  r3 = r2*r;
  A = 0.5*(v[0]*v[0] + v[1]*v[1]) - 1/r;
  eB = -eps*(1.0 - alp*(3*q[0]*q[0])/(r2))/(2*r3);
  return A + eB;
}
