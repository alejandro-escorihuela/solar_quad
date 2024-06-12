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
#define GRAV_CNT 0.000295912208286L

#define IND_Q(I, J, NP) (3*I + J)
#define IND_P(I, J, NP) (3*(I + NP) + J)

real grad(real * z, real * mas, int i, int j, int np);
void kepler_sc(real * z, real * dz, real * mas, int i, real h, int np);
void phiscHK(real * z, real * dz, real * mas, real h, int np);
void phiscHI(real * z, real * dz, real * mas, real h, int np);
real ham(real * z, real * mas, int np);

#endif
