/* -*- mode: C; mode: fold -*- */
/*
 * clib/thresholds.c
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

void thresholds(double* val, int D, 
                const double **arg, const int* len, int N) {
  int Dx, Nth, k, t;
  double s, *z; 
  const double *th, *th_true, *x, *m, *c;
  
  Dx = len[0];
  if (N != 5 || Dx < 1 || len[1] != Dx || len[2] != Dx*Dx || len[3] != D) {
    fprintf(stderr, "thresholds: Invalid arguments.\n");
    exit(EXIT_FAILURE);
  }

  Nth = len[3];
  x = arg[0]; 
  m = arg[1];
  c = arg[2];
  th = arg[3];
  th_true = arg[4];
  
  z = (double *)malloc(sizeof(double)*Dx);
  for (k=0; k<Dx; k++) {
    z[k] = x[k] - m[k];
  }
  tril_solve_inplace(c, z, Dx);
  s = 0.0;
  for (k=0; k<Dx; k++) s += z[k]*z[k];
  for (t=0; t<Nth; t++) {
    val[t] = ((s>th[t]) ? 1.0 : 0.0) - th_true[t];
  }
  free(z);
}
