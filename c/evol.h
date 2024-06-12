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

typedef void (* fluxe) (quad *, quad *, quad *, quad, int);
typedef quad (* qcons) (quad *, quad *, int);

void psc(fluxe phi, quad * z, quad * e, quad * params, quad h, int np);
void pasABA(quad * z, quad * e, quad * params, int np, quad h, quad * a, quad * b, int s, fluxe pA, fluxe pB);
void pasAB(quad * z, quad * e, quad * params, int np, quad h, quad * x, quad * y, int sp, fluxe pA, fluxe pB);
void pasBA(quad * z, quad * e, quad * params, int np, quad h, quad * x, quad * y, int sp, fluxe pA, fluxe pB);

/* Input/Output: */
/*   z: coordenades q i p */
/* Input: */
/*   params: parametres quads, en el cas solar les masses */
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
void evolABAsc(quad * z, quad * params, int np, int Nm, quad h, quad * a, quad * b, int s, quad * x, quad * y, int sp, fluxe pA, fluxe pB, qcons hH, quad * maxerH);
void evolABAsolar_errH(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * lerH);
void evolABAsolar_errQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * lerH);
void evolABAsolar_errHQ(real * z, real * params, int np, int Nm, real h, real * a, real * b, int s, real * x, real * y, int sp, real * lerH, real * lerQ);
#endif
