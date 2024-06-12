/* 10-06-2024 */
/* alex */
/* vector.h */
#ifndef _VECTOR_H
#define _VECTOR_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <quadmath.h>

typedef long double real;
typedef __float128 realq;

real dot(real * v1, real * v2);
void copy(real * v1, real * v2, int tam); /* v1 <- v2 */
void adjunt(real * v1, real * v2, int tam); /* v1 <- adj(v2)*/
real suma(real * v1, int tam);
void zeros(real * v1, int tam);
real difnorm2(real * v1, real * v2, int tam);
void printV(real * v1, const char * nom, int tam);
#endif
