/* 06-06-2024 */
/* alex */
/* solar.c */
#include <stdio.h>
#include <stdlib.h>
#include "solar.h"

quad grad(quad * z, quad * mas, int i, int j, int np) {;
  int k, m;
  quad gV = 0.0, q[np][3], resta[3], den;
  
  if (i == 0)
    return 0.0;
  for (k = 0; k < np; k++)
    for (m = 0; m < 3; m++)
      q[k][m] = z[IND_Q(k, m, np)];
  for (k = 1; k < np; k++)
    if (i != k) {
      for (m = 0; m < 3; m++)
	resta[m] = q[i][m] - q[k][m];
      den = powq((resta[0]*resta[0]) + (resta[1]*resta[1]) + (resta[2]*resta[2]), 1.5);
      gV += (mas[k]*(q[i][j] - q[k][j]))/den;
    }
  gV *= GRAV_CNT*mas[i];
  return gV;  
}

void kepler_sc(quad * z, quad * dz, quad * mas, int i, quad h, int np) {
  /* Sergio Blanes and Fernando Casas: A Concise Introduction to Geometric Numerical Integrator p[28,29] */
  /* Adapted to compensated summation by Ander Murua */
  quad q[3], p[3], q_ant[3], p_ant[3];
  quad t, mu, r0, v02, u, a;
  quad c, s, s2, sig, psi, w, x, dx, adx, adx_ant;
  quad fft, gg, fp, gpt, aux;
  int j;

  if (i > 0) {
    for (j = 0; j < 3; j++) {
      q[j] = q_ant[j] = z[IND_Q(i, j, np)];
      p[j] = p_ant[j] = z[IND_P(i, j, np)];

    }
    t = h/mas[i];
    mu = GRAV_CNT*mas[0]*mas[i]*mas[i];
    r0 = sqrtq(dot(q, q));
    v02 = dot(p, p);
    u = dot(q, p);
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
      dz[IND_Q(i, j, np)] = fft*q_ant[j] + gg*p_ant[j];
      dz[IND_P(i, j, np)] = fp*q_ant[j] + gpt*p_ant[j];
    }
  }
}

void phiscHK(quad * z, quad * dz, quad * mas, quad h, int np) {
  int i;
  
  for (i = 0; i < np; i++)
    kepler_sc(z, dz, mas, i, h, np);
}

void phiscHI(quad * z, quad * dz, quad * mas, quad h, int np) {
  int i, j;
  
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      dz[IND_Q(i, j, np)] = 0.0;
      dz[IND_P(i, j, np)] = -h*grad(z, mas, i, j, np);
    } 
}

quad ham(quad * z, quad * mas, int np) {
  int i, j, k;
  quad cin = 0.0, pot = 0.0, q[np][3], p[np][3], resta[3];
  
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      q[i][j] = z[IND_Q(i, j, np)];
      p[i][j] = z[IND_P(i, j, np)];
    }
  for (i = 0; i < np; i++)
    cin += (p[i][0]*p[i][0] + p[i][1]*p[i][1] + p[i][2]*p[i][2])/mas[i];
  cin *= 0.5;
  for (i = 0; i < np; i++)
    for (j = 0; j < i; j++) {
      for (k = 0; k < 3; k++)
	resta[k] = q[i][k] - q[j][k];
      pot += mas[i]*mas[j]/sqrtq(resta[0]*resta[0] + resta[1]*resta[1] + resta[2]*resta[2]);
    }
  pot *= -GRAV_CNT;
  return cin + pot;
}
