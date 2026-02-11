/* 06-06-2024 */
/* alex */
/* solar.h */
#ifndef _SOLAR_H
#define _SOLAR_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <quadmath.h>
#include "vector.h"

#define IND_Q(I, J, NP) (3*I + J)
#define IND_P(I, J, NP) (3*(I + NP) + J)

void kepler_sc(quad * z, quad * dz, int nz, quad * par, int i, quad h, int np);
void phiscHK(quad * z, quad * dz, int nz, quad * par, quad h, int np);
void phiscHI(quad * z, quad * dz, int nz, quad * par, quad h, int np);
quad ham(quad * z, int nz, quad * par, int np);

void expand_masses(quad * par, quad * GM, int np);
void centrar(quad * z, quad * par, int np);
void cart2jacobi(quad * zb, quad * z, quad * par, int np);
void jacobi2cart(quad * z, quad * zb, quad * par, int np);
#endif
