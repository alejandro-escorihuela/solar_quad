/* 06-06-2024 */
/* alex */
/* solar.c */
#include <stdio.h>
#include <stdlib.h>
#include "solar.h"

void kepler_sc(quad * z, quad * dz, quad * par, int i, quad h, int np) {
  /* Sergio Blanes and Fernando Casas: A Concise Introduction to Geometric Numerical Integrator p[28,29] */
  /* Adapted to compensated summation by Ander Murua */
  quad q[3], v[3], q_ant[3], v_ant[3];
  quad t, mu, r0, v02, u, a;
  quad c, s, s2, sig, psi, w, x, dx, adx, adx_ant;
  quad fft, gg, fp, gpt, aux;
  int j;

  if (i > 0) {
    for (j = 0; j < 3; j++) {
      q[j] = q_ant[j] = z[IND_Q(i, j, np)];
      v[j] = v_ant[j] = z[IND_P(i, j, np)];

    }
    t = h;
    mu = par[np + i];
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
      dz[IND_Q(i, j, np)] = fft*q_ant[j] + gg*v_ant[j];
      dz[IND_P(i, j, np)] = fp*q_ant[j] + gpt*v_ant[j];
    }
  }
}

void phiscHK(quad * z, quad * dz, quad * par, quad h, int np) {
  int i;
  
  for (i = 1; i < np; i++)
    kepler_sc(z, dz, par, i, h, np);
}

void phiscHI(quad * z, quad * dz, quad * par, quad h, int np) {
  int i, j, k;
  quad nu[np], mu[np];  
  quad q[np][3], Q[np][3], qpp[np][3], qbpp[np][3];
  quad resta[3], mod, den;
  
  /* par -> m, nu */
  for (i = 0; i < np; i++) {
    mu[i] = par[np + i];
    nu[i] = par[2*np + i];
  }
  
  /* qb -> q */
  for (j = 0; j < 3; j++)
    Q[np - 1][j] = z[IND_Q(0, j, np)];
  for (i = np - 1; i > 0; i--)
    for (j = 0; j < 3; j++) {
      Q[i - 1][j] = Q[i][j] - nu[i]*z[IND_Q(i, j, np)];
      q[i][j] = z[IND_Q(i, j, np)] + Q[i - 1][j];   
    }
  for (j = 0; j < 3; j++)
    q[0][j] = Q[0][j];
  /* qpp */
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++)
      qpp[i][j] = 0.0;
  for (i = 0; i < np; i++)
    for (k = 0; k < np; k++)
      if (i != k) {
	for (j = 0; j < 3; j++)
	  resta[j] = q[i][j] - q[k][j];
	den = 1.0/powq((resta[0]*resta[0]) + (resta[1]*resta[1]) + (resta[2]*resta[2]), 1.5);
	for (j = 0; j < 3; j++)
	  qpp[i][j] -= par[k]*resta[j]*den;
      }
  /* qpp -> qbpp */
  for (j = 0; j < 3; j++)
    Q[0][j] = qpp[0][j];
  for (i = 1; i < np; i++)
    for (j = 0; j < 3; j++) {
      qbpp[i][j] = qpp[i][j] - Q[i - 1][j];
      Q[i][j] = Q[i - 1][j] + nu[i]*qbpp[i][j];  
    }
  for (j = 0; j < 3; j++)
    qbpp[0][j] = Q[np - 1][j];
  /* dvb */
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++)
      dz[IND_Q(i, j, np)] = dz[IND_P(i, j, np)] = 0.0;
  for (i = 1; i < np; i++) {
    mod = 1.0/powq(z[IND_Q(i, 0, np)]*z[IND_Q(i, 0, np)] + z[IND_Q(i, 1, np)]*z[IND_Q(i, 1, np)] + z[IND_Q(i, 2, np)]*z[IND_Q(i, 2, np)], 1.5);
    for (j = 0; j < 3; j++)
      dz[IND_P(i, j, np)] = h*(mu[i]*z[IND_Q(i, j, np)]*mod + qbpp[i][j]);
  }
}

quad ham(quad * z, quad * par, int np) {
  int i, j, k;
  quad zcart[3*np*2], cin = 0.0, pot = 0.0, q[np][3], v[np][3], resta[3];
  
  jacobi2cart(zcart, z, par, np);
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      q[i][j] = zcart[IND_Q(i, j, np)];
      v[i][j] = zcart[IND_P(i, j, np)];
    }
  for (i = 0; i < np; i++)
    cin += par[i]*(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  cin *= 0.5;
  for (i = 1; i < np; i++)
    for (j = 0; j < i; j++) {
      for (k = 0; k < 3; k++)
	resta[k] = q[i][k] - q[j][k];
      pot += par[i]*par[j]/sqrtq(resta[0]*resta[0] + resta[1]*resta[1] + resta[2]*resta[2]);
    }
  return cin - pot;
}

void expand_masses(quad * par, quad * Gm, int np) {
  int i;
  quad GM[np], nu[np], mu[np];
  
  GM[0] = Gm[0];
  for (i = 1; i < np; i++)
    GM[i] = GM[i - 1] + Gm[i];
  nu[0] = Gm[0]/GM[0];
  mu[0] = 0.0;
  for (i = 1; i < np; i++) {
    mu[i] = Gm[0]*GM[i]/GM[i - 1];
    nu[i] = Gm[i]/GM[i];
  }
  for (i = 0; i < np; i++) {
    par[i] = Gm[i];
    par[np + i] = mu[i];
    par[2*np + i] = nu[i];
  }
}

void centrar(quad * z, quad * par, int np) {
  int i, j;
  quad M = 0, qcm[3] = {0, 0, 0}, vcm[3] = {0, 0, 0};

  for (i = 0; i < np; i++) {
    M += par[i];
    for (j = 0; j < 3; j++) {
      qcm[j] += par[i]*z[IND_Q(i, j, np)];
      vcm[j] += par[i]*z[IND_P(i, j, np)];
    }
  }
  for (j = 0; j < 3; j++) {
    qcm[j] /= M;
    vcm[j] /= M;
  }
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      z[IND_Q(i, j, np)] -= qcm[j];
      z[IND_P(i, j, np)] -= vcm[j];
    }  
}

void cart2jacobi(quad * zb, quad * z, quad * par, int np) {
  int i, j;
  quad nu[np];
  quad qb[np][3], Q[np][3];
  quad vb[np][3], V[np][3];
  
  /* par -> m, nu */
  for (i = 0; i < np; i++)
    nu[i] = par[2*np + i];    

  /* q, v -> qb, vb */
  for (j = 0; j < 3; j++) {
    Q[0][j] = z[IND_Q(0, j, np)];
    V[0][j] = z[IND_P(0, j, np)];
  }
  for (i = 1; i < np; i++)
    for (j = 0; j < 3; j++) {
      qb[i][j] = z[IND_Q(i, j, np)] - Q[i - 1][j];
      Q[i][j] = Q[i - 1][j] + nu[i]*qb[i][j];
      vb[i][j] = z[IND_P(i, j, np)] - V[i - 1][j]; 
      V[i][j] = V[i - 1][j] + nu[i]*vb[i][j];      
    }
  for (j = 0; j < 3; j++) {
    qb[0][j] = Q[np - 1][j];
    vb[0][j] = V[np - 1][j];      
  }

  /* qb, vb -> zb */
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      zb[IND_Q(i, j, np)] = qb[i][j];
      zb[IND_P(i, j, np)] = vb[i][j];
    }
}

void jacobi2cart(quad * z, quad * zb, quad * par, int np) {
  int i, j;
  quad nu[np];
  quad q[np][3], Q[np][3];
  quad v[np][3], V[np][3];
  
  /* par -> m, nu */
  for (i = 0; i < np; i++)
    nu[i] = par[2*np + i];    
  
  /* qb, vb -> q, v */
  for (j = 0; j < 3; j++) {
    Q[np - 1][j] = zb[IND_Q(0, j, np)];
    V[np - 1][j] = zb[IND_P(0, j, np)];
  }
  for (i = np - 1; i > 0; i--)
    for (j = 0; j < 3; j++) {
      Q[i - 1][j] = Q[i][j] - nu[i]*zb[IND_Q(i, j, np)];
      q[i][j] = zb[IND_Q(i, j, np)] + Q[i - 1][j];
      V[i - 1][j] = V[i][j] - nu[i]*zb[IND_P(i, j, np)];
      v[i][j] = zb[IND_P(i, j, np)] + V[i - 1][j]; 
    }
  for (j = 0; j < 3; j++) {
    q[0][j] = Q[0][j];
    v[0][j] = V[0][j];
  }

  /* q, v -> z */
  for (i = 0; i < np; i++)
    for (j = 0; j < 3; j++) {
      z[IND_Q(i, j, np)] = q[i][j];
      z[IND_P(i, j, np)] = v[i][j];
    }
}
