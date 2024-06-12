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
#undef GRAV_CNT
#define GRAV_CNT 0.000295912208286q

#define IND_Q(I, J, NP) (3*I + J)
#define IND_P(I, J, NP) (3*(I + NP) + J)

quad grad(quad * z, quad * mas, int i, int j, int np);
void kepler_sc(quad * z, quad * dz, quad * mas, int i, quad h, int np);
void phiscHK(quad * z, quad * dz, quad * mas, quad h, int np);
void phiscHI(quad * z, quad * dz, quad * mas, quad h, int np);
quad ham(quad * z, quad * mas, int np);

#endif
