/* -*- mode: C; mode: fold -*- */
/*
 * clib/dbanana.c
 * 
 * Copyright (c) Matti Vihola 2009-2013
 * 
 * This file is part of Grapham.
 * 
 * Grapham is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Grapham is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Grapham.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include "grapham_math.h"

double dbanana(const double **arg, const int* len, const int N) {
  double tmp[2];
  double a, b;
  const double *x;
  if (N != 5 || len[0] != 2) {
    fprintf(stderr, "dbanana: Invalid arguments.\n");
    exit(EXIT_FAILURE);
  }
  x = arg[0];
  a = arg[3][0];
  b = arg[4][0];
  tmp[0] = x[0]/a;
  tmp[1] = x[1]*a + a*b*(x[0]*x[0] + a*a);
  return(d_mvnorm_chol(tmp, arg[1], arg[2], 2));
}

void banana_fun(double* val, int D, 
                const double **arg, const int* len, int N) {
  const double pp[] = {0.21072103131565,
   0.57536414490356,
   1.38629436111989,
   2.77258872223979,
   4.60517018598809};
  const double true_pp[] = {.1, .25, .5, .75, .9};
  const int Npp = 5;
  
  double tmp[2];
  double z, a, b;
  int k;
  const double *x, *v;
  
  x = arg[0];
  v = arg[1];
  a = arg[2][0];
  b = arg[3][0];
  
  tmp[0] = x[0]/a;
  tmp[1] = x[1]*a + a*b*(x[0]*x[0] + a*a);
  
  tril_solve_inplace(v, tmp, 2);
  
  z = 0;
  for (k=0; k<2; k++) z += tmp[k]*tmp[k];
  
  for (k=0; k<Npp; k++) val[k] = 1-true_pp[k];
  for (k=0; k<Npp; k++) {
     if (z<=pp[k]) {
       break;
     } else {
       val[k] = val[k]-1;
     }
  }
}
