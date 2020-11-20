/* -*- mode: C; mode: fold -*- */
/*
 * clib/dcircular.c
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

double dcircular(const double **arg, const int* len, const int N) {
  if (N != 2 || len[0] != 2 || len[1] != 1) {
    fprintf(stderr, "dcircular: Invalid arguments.\n");
    exit(EXIT_FAILURE);
  }
  const double* x = arg[0];
  double c = arg[1][0];
  double nsqx = c-sqrt(x[0]*x[0]+x[1]*x[1]);
  return -nsqx*nsqx + log((1+cos(5*x[0]))*(1+cos(5*x[1])));
}
