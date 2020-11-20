/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_lib.h
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

#ifndef _GRAPHAM_LIB_H
#define _GRAPHAM_LIB_H

typedef struct {
  char *name;
  double (*c_fun)(double **, int *, int);
  int (*L_fun)(lua_State *);
} grapham_dlib_struct;

typedef struct {
  char *name;
  void (*c_fun)(double*, int, const double **, const int *, int);
} grapham_flib_struct;

#ifndef _GRAPHAM_LIB_C
extern grapham_dlib_struct grapham_dlib[];
extern grapham_flib_struct grapham_flib[];
#endif

#define GRAPHAM_LIB_DISTPROTO(dist) \
  double d##dist(double **, int *, int)

GRAPHAM_LIB_DISTPROTO(beta);
GRAPHAM_LIB_DISTPROTO(binom); 
GRAPHAM_LIB_DISTPROTO(cauchy);
GRAPHAM_LIB_DISTPROTO(chi2);
GRAPHAM_LIB_DISTPROTO(erlang);
GRAPHAM_LIB_DISTPROTO(exp); 
GRAPHAM_LIB_DISTPROTO(fisher); 
GRAPHAM_LIB_DISTPROTO(gamma); 
GRAPHAM_LIB_DISTPROTO(gumbel);
GRAPHAM_LIB_DISTPROTO(invchi2); 
GRAPHAM_LIB_DISTPROTO(invgamma); 
GRAPHAM_LIB_DISTPROTO(laplace); 
GRAPHAM_LIB_DISTPROTO(levy);
GRAPHAM_LIB_DISTPROTO(logistic); 
GRAPHAM_LIB_DISTPROTO(mvnorm);
GRAPHAM_LIB_DISTPROTO(mvnorm_chol);
GRAPHAM_LIB_DISTPROTO(mvnorm_diag);
GRAPHAM_LIB_DISTPROTO(mvstudent);
GRAPHAM_LIB_DISTPROTO(nbinom); 
GRAPHAM_LIB_DISTPROTO(norm); 
GRAPHAM_LIB_DISTPROTO(pareto); 
GRAPHAM_LIB_DISTPROTO(poisson);
GRAPHAM_LIB_DISTPROTO(rayleigh);
GRAPHAM_LIB_DISTPROTO(student); 
GRAPHAM_LIB_DISTPROTO(uniform);
GRAPHAM_LIB_DISTPROTO(wishart); 

#define GRAPHAM_LIB_FUNCPROTO(fun) \
  void f##fun(double*, int, const double **, const int *, int)

GRAPHAM_LIB_FUNCPROTO(mean);
GRAPHAM_LIB_FUNCPROTO(mom2);

#endif
