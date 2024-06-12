/* 06-06-2024 */
/* alex */
/* solar.c */
#include <stdio.h>
#include <stdlib.h>
#include "solar.h"

real grad(real * z, real * mas, int i, int j, int np) {;
  int k, m;
  real gV = 0.0, q[np][3], resta[3], den;
  
  if (i == 0)
    return 0.0L;
  for (k = 0; k < np; k++)
    for (m = 0; m < 3; m++)
      q[k][m] = z[IND_Q(k, m, np)];
  for (k = 1; k < np; k++)
    if (i != k) {
      for (m = 0; m < 3; m++)
	resta[m] = q[i][m] - q[k][m];
      den = powl((resta[0]*resta[0]) + (resta[1]*resta[1]) + (resta[2]*resta[2]), 1.5L);
      gV += (mas[k]*(q[i][j] - q[k][j]))/den;
    }
  gV *= GRAV_CNT*mas[i];
  return gV;  
}

void kepler_sc(real * z, real * dz, real * mas, int i, real h, int np) {
  /* Sergio Blanes and Fernando Casas: A Concise Introduction to Geometric Numerical Integrator p[28,29] */
  /* Adapted to compensated summation by Ander Murua */
  real q[3], p[3], q_ant[3], p_ant[3];
  real t, mu, r0, v02, u, a;
  real c, s, s2, sig, psi, w, x, dx, adx, adx_ant;
  real fft, gg, fp, gpt, aux;
  int j;

  if (i > 0) {
    for (j = 0; j < 3; j++) {
      q[j] = q_ant[j] = z[IND_Q(i, j, np)];
      p[j] = p_ant[j] = z[IND_P(i, j, np)];

    }
    t = h/mas[i];
    mu = GRAV_CNT*mas[0]*mas[i]*mas[i];
    r0 = sqrtl(dot(q, q));
    v02 = dot(p, p);
    u = dot(q, p);
    a = -mu/(v02 - ((2.0L*mu)/r0));
    w = sqrtl(mu/(a*a*a));
    sig = 1.0L - r0/a;
    psi = u/(w*a*a);
    adx = INFINITY;
    x = w*t*(a/r0);
    do {
      adx_ant = adx;
      c = cosl(x);
      s = sinl(x);
      s2 = sinl(x/2.0L);
      s2 = 2.0L*s2*s2;
      dx = (x - sig*s + psi*s2 - w*t)/(1.0L - sig*c + psi*s);
      x = x - dx;
      adx = fabsl(dx);
    } while (adx != 0.0L && adx < adx_ant);
    aux = 1.0L - sig*c + psi*s;
    fft = -s2*a/r0;
    gg = t + ((s - x)/w);
    fp = (-a*w*s)/(aux*r0);
    gpt = -s2/aux;
    for (j = 0; j < 3; j++) {
      dz[IND_Q(i, j, np)] = fft*q_ant[j] + gg*p_ant[j];
      dz[IND_P(i, j, np)] = fp*q_ant[j] + gpt*p_ant[j];
    }
  }
}

void phiscHK(real * z, real * dz, real * mas, real h, int np) {
  int i;
  
  for (i = 0; i < np; i++)
    kepler_sc(z, dz, mas, i, h, np);
}

void phiscHI(real * z, real * dz, real * mas, real h, int np) {
  int i, j;
  
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      dz[IND_Q(i, j, np)] = 0.0L;
      dz[IND_P(i, j, np)] = -h*grad(z, mas, i, j, np);
    } 
}

real ham(real * z, real * mas, int np) {
  int i, j, k;
  real cin = 0.0L, pot = 0.0L, q[np][3], p[np][3], resta[3];
  
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      q[i][j] = z[IND_Q(i, j, np)];
      p[i][j] = z[IND_P(i, j, np)];
    }
  for (i = 0; i < np; i++)
    cin += (p[i][0]*p[i][0] + p[i][1]*p[i][1] + p[i][2]*p[i][2])/mas[i];
  cin *= 0.5L;
  for (i = 0; i < np; i++)
    for (j = 0; j < i; j++) {
      for (k = 0; k < 3; k++)
	resta[k] = q[i][k] - q[j][k];
      pot += mas[i]*mas[j]/sqrtl(resta[0]*resta[0] + resta[1]*resta[1] + resta[2]*resta[2]);
    }
  pot *= -GRAV_CNT;
  return cin + pot;
}
