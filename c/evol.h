/* 10-06-2024 */
/* alex */
/* evol.h */
#ifndef _EVOL_H
#define _EVOL_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <quadmath.h>
#include "vector.h"
#include "solar.h"

typedef void (* fluxe) (real *, real *, real *, real, int);
typedef real (* qcons) (real *, real *, int);

void psc(fluxe phi, real * z, real * e, real * params, real h, int np);
void pasABA(real * z, real * e, real * params, int np, real h, real * a, real * b, int s, fluxe pA, fluxe pB);
void pasAB(real * z, real * e, real * params, int np, real h, real * x, real * y, int sp, fluxe pA, fluxe pB);
void pasBA(real * z, real * e, real * params, int np, real h, real * x, real * y, int sp, fluxe pA, fluxe pB);
void evolABAsc(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, fluxe pA, fluxe pB, qcons hH, real * maxerH);
/* Input/Output: */
/*   z: coordenades q i p */
/* Input: */
/*   params: parametres reals, en el cas solar les masses */
/*   np: dimensió de params, en el cas solar nombre de planetes */
/*   h: grandària de pas */
/*   a, b: coeficients a i b */
/*   x, y: coeficients a i b del processat */
/*   s, sp: nombre d'etapes del mètode i del processat */
/*   pA, pB: fluxes A i B amb suma compensada */
/*   hH: funció de la quantitat conservada */
/* Output: */
/*   maxerH: màxim error de la quatitat conservada */
/* Queda pendent escriure en fitxers per veure l'evolució de les z */

void evolABAsolar_errH(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * erH);
void evolABAsolar_errQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * erH);
void evolABAsolar_errHQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * erH, real * erQ);
#endif
