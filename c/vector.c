/* 10-06-2024 */
/* alex */
/* vector.c */
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

quad dot(quad * v1, quad * v2) {
  int i;
  quad ret = 0.0;
  for (i = 0; i < 3; i++)
    ret += v1[i]*v2[i];
  return ret;
}

void copy(quad * v1, quad * v2, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = v2[i];
}

void adjunt(quad * v1, quad * v2, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = -v2[tam - 1 - i];
}

quad suma(quad * v1, int tam) {
  int i;
  quad s = 0.0;
  for (i = 0; i < tam; i++)
    s += v1[i];
  return s;
}

void zeros(quad * v1, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = 0.0;
}

quad difnorm2(quad * v1, quad * v2, int tam) {
  int i;
  quad nu = 0.0, de = 0.0;
  for (i = 0; i < tam; i++) {
    nu += (v1[i] - v2[i])*(v1[i] - v2[i]);
    de += v2[i]*v2[i];
  }
  return sqrtq(nu/de);
}

void printV(quad * v1, const char * nom, int tam) {
  int i;
  printf("%s = {", nom);
  for (i = 0; i < tam - 1; i++) {
    if (i % 5 == 0)
      printf("\n  ");    
    printf("%02d => %Qe, ", i, v1[i]);

  }
  if ((tam - 1) % 5 == 0)
    printf("\n  "); 
  printf("%02d => %Qe\n}\n", tam - 1, v1[tam - 1]);
  
}

quad real2quad(real v1) {
  char num[26];
  sprintf(num, "%.19Le", v1);
  //printf("%s\n", num);
  return strtoflt128 (num, NULL);
}

real quad2real(quad v1) {
  char num[41];
  sprintf(num, "%.34Qe", v1);
  return atof(num);
}

void real2quadV(quad * v1, real * v2, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = real2quad(v2[i]);  
}

void quad2realV(real * v1, quad * v2, int tam) {
  int i;
  for (i = 0; i < tam; i++)
    v1[i] = quad2real(v2[i]);
}
