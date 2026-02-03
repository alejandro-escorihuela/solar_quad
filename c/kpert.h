/* 03-02-2026 */
/* alex */
/* kpert.h */
#ifndef _KPERT_H
#define _KPERT_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <quadmath.h>
#include "vector.h"

void kepler_sc_kp(quad * z, quad * dz, quad h);
void phiscH0_kp(quad * z, quad * dz, quad * par, quad h, int np);
void phiscH1_kp(quad * z, quad * dz, quad * par, quad h, int np);
quad ham_kp(quad * z, quad * par, int np);

#endif
