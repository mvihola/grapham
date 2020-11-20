/* -*- mode: C; mode: fold -*- */
/*
 * clib/dcross.c
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

double dcross(const double **arg, const int* len, const int N) {
  double p1, p2;
  const double 
  m1[] = {4,4};
  const double
  m2[] = {0,0},
  CC1[] = {1.00249688278817, 0,
           -0.99252178942710, 0.14106912317172},
  CC2[] = {1.00249688278817, 0,
           0.99252178942710, 0.14106912317172};
  
  if (N != 1 || len[0] != 2) {
    fprintf(stderr, "dcircular: Invalid arguments.\n");
    exit(EXIT_FAILURE);
  }
  
  p1 = d_mvnorm_chol(arg[0], m1, CC1, 2);
  p2 = d_mvnorm_chol(arg[0], m2, CC2, 2);
  
  return(log(exp(p1)+exp(p2)));
}

void cross_fun(double* val, int D, 
                const double **arg, const int* len, int N) {
  if (D != 2 || N != 1 || len[0] != 2) {
    fprintf(stderr, "cross_fun: Invalid arguments.\n");
    exit(EXIT_FAILURE);
  }
  
  val[0] = arg[0][0] - 2;
  val[1] = arg[0][1] - 2;
}
