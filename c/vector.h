/* 10-06-2024 */
/* alex */
/* vector.h */
#ifndef _VECTOR_H
#define _VECTOR_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <quadmath.h>

typedef long double real;
typedef __float128 quad;

quad dot(quad * v1, quad * v2, int tam);
void copy(quad * v1, quad * v2, int tam);   /* v1 <- v2 */
void adjunt(quad * v1, quad * v2, int tam); /* v1 <- adj(v2)*/
quad suma(quad * v1, int tam);
void zeros(quad * v1, int tam);
quad difnorm2(quad * v1, quad * v2, int tam);
void printV(quad * v1, const char * nom, int tam);
void printV_str(const char ** v1, const char * nom, int tam);

quad real2quad(real v1);
real quad2real(quad v1);
void real2quadV(quad * v1, real * v2, int tam); /* v1 <- quad(v2) */
void quad2realV(real * v1, quad * v2, int tam); /* v1 <- real(v2) */
quad str2quad(const char * v1);
const char * quad2str(quad v1);
void str2quadV(quad * v1, const char ** v2, int tam);
void quad2strV(char ** v1, quad * v2, int tam);
#endif
