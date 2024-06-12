/* 10-06-2024 */
/* alex */
/* vector.c */
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

real dot(real * v1, real * v2) {
  int i;
  real ret = 0.0;
  for (i = 0; i < 3; i++)
    ret += v1[i]*v2[i];
  return ret;
}

void copy(real * v1, real * v2, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = v2[i];
}

void adjunt(real * v1, real * v2, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = -v2[tam - 1 - i];
}

real suma(real * v1, int tam) {
  int i;
  real s = 0.0L;
  for (i = 0; i < tam; i++)
    s += v1[i];
  return s;
}

void zeros(real * v1, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = 0.0L;
}

real difnorm2(real * v1, real * v2, int tam) {
  int i;
  real nu = 0.0L, de = 0.0L;
  for (i = 0; i < tam; i++) {
    nu += (v1[i] - v2[i])*(v1[i] - v2[i]);
    de += v2[i]*v2[i];
  }
  return sqrtl(nu/de);
}

void printV(real * v1, const char * nom, int tam) {
  int i;
  printf("%s = {", nom);
  for (i = 0; i < tam - 1; i++) {
    if (i % 5 == 0)
      printf("\n  ");    
    printf("%02d => %Le, ", i, v1[i]);

  }
  if ((tam - 1) % 5 == 0)
    printf("\n  "); 
  printf("%02d => %Le\n}\n", tam - 1, v1[tam - 1]);
  
}
