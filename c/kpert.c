/* 03-02-2026 */
/* alex */
/* kpert.c */
#include <stdio.h>
#include <stdlib.h>
#include "kpert.h"

void kepler_sc_kp(quad * z, quad * dz, quad h) {
  /* Sergio Blanes and Fernando Casas: A Concise Introduction to Geometric Numerical Integrator p[28,29] */
  /* Adapted to compensated summation by Ander Murua */
  quad q[3], v[3];
  quad t, mu = 1.0Q, r0, v02, u, a;
  quad c, s, s2, sig, psi, w, x, dx, adx, adx_ant;
  quad fft, gg, fp, gpt, aux;
  int j;

  for (j = 0; j < 3; j++) {
    q[j] = z[j];
    v[j] = z[j + 3];  
  }
  t = h;
  r0 = sqrtq(dot(q, q));
  v02 = dot(v, v);
  u = dot(q, v);
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
  for (j = 0; j < 3; j++) {
    dz[j] = fft*q[j] + gg*v[j];
    dz[j + 3] = fp*q[j] + gpt*v[j];
  }
}

void phiscH0_kp(quad * z, quad * dz, quad * par, quad h, int np) {
  (void) par;
  (void) np;
  kepler_sc_kp(z, dz, h);
}
void phiscH1_kp(quad * z, quad * dz, quad * par, quad h, int np);
quad ham_kp(quad * z, quad * par, int np);
