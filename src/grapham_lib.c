/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_lib.c
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

#define _GRAPHAM_LIB_C
#include <string.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#ifdef _HAVE_NUMLUA
#include <luamatrix.h>
#endif

#include "lua_tools.h"
#include "grapham_math.h"
#include "grapham_lib.h"

/* To avoid redundancy, the C and Lua wrapper function definitions 
 * are macros */

/* Simple two and three argument functions */ /*{{{*/
#define D_TWOARG(dist) \
  double d##dist(double **arg, int *len, int N) { \
  if (N != 2 || len[0] != 1 || len[1] != 1) { \
    fprintf(stderr, "d" #dist ": Invalid arguments.\n"); \
    exit(EXIT_FAILURE); \
  }  \
  return d_##dist(arg[0][0], arg[1][0]); \
} \
static int L_d##dist(lua_State* L) { \
  double a0, a1; \
  a0 = luaL_checknumber(L, 1); \
  a1 = luaL_checknumber(L, 2); \
  lua_pushnumber(L, d_##dist(a0, a1)); \
  return 1; \
}


#define D_THREEARG(dist) \
  double d##dist(double **arg, int *len, int N) { \
  if (N != 3 || len[0] != 1 || len[1] != 1 || len[2] != 1) { \
    fprintf(stderr, "d" #dist ": Invalid arguments.\n"); \
    exit(EXIT_FAILURE); \
  }  \
  return d_##dist(arg[0][0], arg[1][0], arg[2][0]); \
} \
static int L_d##dist(lua_State* L) { \
  double a0, a1, a2; \
  a0 = luaL_checknumber(L, 1); \
  a1 = luaL_checknumber(L, 2); \
  a2 = luaL_checknumber(L, 3); \
  lua_pushnumber(L, d_##dist(a0, a1, a2)); \
  return 1; \
}
/*}}}*/

/* Multidimensional */ /*{{{*/
#define D_VECVECMAT(dist) \
double d##dist(double **arg, int* len, int N) { \
  int D; \
  double res, *a0, *a1, *a2; \
  D = len[0]; \
  if (N != 3 || len[1]!=D || len[2]!=D*D) { \
    fprintf(stderr, "d" #dist ": Invalid arguments.\n"); \
    exit(EXIT_FAILURE); \
  } \
  a0 = (double *)malloc(sizeof(double)*(2*D+D*D)); \
  a1 = &(a0[D]); a2 = &(a0[2*D]); \
  memcpy(a0, arg[0], sizeof(double)*D); \
  memcpy(a1, arg[1], sizeof(double)*D); \
  memcpy(a2, arg[2], sizeof(double)*D*D); \
  res = d_##dist(a0,a1,a2,D); \
  free(a0); \
  return res; \
} \
static int L_d##dist(lua_State* L) { \
  double p, *a0, *a1, *a2; \
  int D; \
  D = lua_getlength(L, 1); \
  if (D<0) error(L, "Error: d" #dist ": invalid arguments"); \
  a0 = (double *)malloc(sizeof(double)*(2*D + D*D)); \
  a1 = &(a0[D]); a2 = &(a0[2*D]); \
  lua_getvector(L, 1, a0, D); \
  lua_getvector(L, 2, a1, D); \
  lua_gettril(L, 3, a2, D); \
  p = d_##dist(a0, a1, a2, D); \
  free(a0); \
  lua_pushnumber(L, p); \
  return 1; \
}

#define D_VECVECVEC(dist) \
double d##dist(double **arg, int* len, int N) { \
  int D; \
  double res, *a0, *a1, *a2; \
  D = len[0]; \
  if (N != 3 || len[1]!=D || len[2]!=D) { \
    fprintf(stderr, "d" #dist ": Invalid arguments.\n"); \
    exit(EXIT_FAILURE); \
  } \
  a0 = (double *)malloc(sizeof(double)*(3*D)); \
  a1 = &(a0[D]); a2 = &(a0[2*D]); \
  memcpy(a0, arg[0], sizeof(double)*D); \
  memcpy(a1, arg[1], sizeof(double)*D); \
  memcpy(a2, arg[2], sizeof(double)*D); \
  res = d_##dist(a0,a1,a2,D); \
  free(a0); \
  return res; \
} \
static int L_d##dist(lua_State* L) { \
  double p, *a0, *a1, *a2; \
  int D; \
  D = lua_getlength(L, 1); \
  if (D<0) error(L, "Error: d" #dist ": invalid arguments"); \
  a0 = (double *)malloc(sizeof(double)*(3*D)); \
  a1 = &(a0[D]); a2 = &(a0[2*D]); \
  lua_getvector(L, 1, a0, D); \
  lua_getvector(L, 2, a1, D); \
  lua_getvector(L, 3, a2, D); \
  p = d_##dist(a0, a1, a2, D); \
  free(a0); \
  lua_pushnumber(L, p); \
  return 1; \
}

#define D_VECVECMATSC(dist) \
double d##dist(double **arg, int* len, int N) { \
  int D; \
  double res, *a0, *a1, *a2; \
  D = len[0]; \
  if (N != 4 || len[1]!=D || len[2]!=D*D || len[3]!=1) { \
    fprintf(stderr, "d" #dist ": Invalid arguments.\n"); \
    exit(EXIT_FAILURE); \
  } \
  a0 = (double *)malloc(sizeof(double)*(2*D+D*D)); \
  a1 = &(a0[D]); a2 = &(a0[2*D]); \
  memcpy(a0, arg[0], sizeof(double)*D); \
  memcpy(a1, arg[1], sizeof(double)*D); \
  memcpy(a2, arg[2], sizeof(double)*D*D); \
  res = d_##dist(a0,a1,a2,arg[3][0],D); \
  free(a0); \
  return res; \
} \
static int L_d##dist(lua_State* L) { \
  double p, *a0, *a1, *a2, a3; \
  int D; \
  D = lua_getlength(L, 1); \
  if (D<0) error(L, "Error: d" #dist ": invalid arguments"); \
  a0 = (double *)malloc(sizeof(double)*(2*D + D*D)); \
  a1 = &(a0[D]); a2 = &(a0[2*D]); \
  lua_getvector(L, 1, a0, D); \
  lua_getvector(L, 2, a1, D); \
  lua_gettril(L, 3, a2, D); \
  a3 = luaL_checknumber(L, 4); \
  p = d_##dist(a0, a1, a2, a3, D); \
  free(a0); \
  lua_pushnumber(L, p); \
  return 1; \
}

/*}}}*/

/* Continuous, simple */
D_THREEARG(beta)
D_THREEARG(cauchy)
D_TWOARG(chi2)
D_THREEARG(erlang)
D_TWOARG(exp)
D_THREEARG(fisher)
D_THREEARG(gamma)
D_THREEARG(gumbel)
D_TWOARG(invchi2)
D_THREEARG(invgamma)
D_THREEARG(laplace)
D_TWOARG(levy)
D_THREEARG(logistic)
D_THREEARG(lognorm)
D_THREEARG(norm)
D_THREEARG(pareto)
D_TWOARG(rayleigh)
D_TWOARG(student)

/* Discrete */
D_THREEARG(binom)
D_THREEARG(nbinom)
D_TWOARG(poisson)

/* Multidimensional */
D_VECVECMAT(mvnorm)
D_VECVECMAT(mvnorm_chol)
D_VECVECVEC(mvnorm_diag)
D_VECVECMATSC(mvstudent)

/* Dummy "uniform distribution" returning zero independently how
 * called. */
double duniform(double **arg, int *len, int N) { /*{{{*/
  return 0.0;
} /*}}}*/
static int L_duniform(lua_State* L) { /*{{{*/
  lua_pushnumber(L, 0.0);
  return 1;
} /*}}}*/

double dwishart(double **arg, int *len, int N) { /*{{{*/
  int D, L;
  double res, Dsqrt, *x, *v;
  
  L = len[0];
  Dsqrt = sqrt(L);
  D = Dsqrt;
  if (N != 3 || (double)D != Dsqrt || len[1]!=L || len[2]!=1) { 
    fprintf(stderr, "dwishart: Invalid arguments.\n"); 
    exit(EXIT_FAILURE); 
  } 
  x = (double *)malloc(sizeof(double)*(2*L)); 
  v = &(x[L]); 
  memcpy(x, arg[0], sizeof(double)*L); 
  memcpy(v, arg[1], sizeof(double)*L); 
  res = d_wishart(x,v,arg[2][0],D); 
  free(x); 
  return res; 
} /*}}}*/
static int L_dwishart(lua_State* L) { /*{{{*/
  double p, n;
  double *x, *v;
  int N;
  
  N = lua_getlength(L, 1);
  if (N<0) error(L, "Error: d_wishart: invalid arguments");
  
  x = (double *)malloc(sizeof(double)*(2*N*N));
  v = &(x[N*N]); 
  
  lua_gettril(L, 1, x, N);
  lua_gettril(L, 2, v, N);
  n = luaL_checknumber(L, 3);
  
  p = d_wishart(x, v, n, N);
  
  free(x); 
  
  lua_pushnumber(L, p);
  return 1;
}

/*}}}*/

/* Define a look up table for both C and Lua functions. */
#define GRAPHAM_DLIB_LUT(dist) {"d" #dist, &d##dist, &L_d##dist}

const grapham_dlib_struct grapham_dlib[] = {
  GRAPHAM_DLIB_LUT(beta),
  GRAPHAM_DLIB_LUT(binom), 
  GRAPHAM_DLIB_LUT(cauchy),
  GRAPHAM_DLIB_LUT(chi2),
  GRAPHAM_DLIB_LUT(erlang),
  GRAPHAM_DLIB_LUT(exp), 
  GRAPHAM_DLIB_LUT(fisher), 
  GRAPHAM_DLIB_LUT(gamma), 
  GRAPHAM_DLIB_LUT(gumbel),
  GRAPHAM_DLIB_LUT(invchi2), 
  GRAPHAM_DLIB_LUT(invgamma), 
  GRAPHAM_DLIB_LUT(laplace), 
  GRAPHAM_DLIB_LUT(levy),
  GRAPHAM_DLIB_LUT(logistic), 
  GRAPHAM_DLIB_LUT(lognorm), 
  GRAPHAM_DLIB_LUT(mvnorm),
  GRAPHAM_DLIB_LUT(mvnorm_chol),
  GRAPHAM_DLIB_LUT(mvnorm_diag),
  GRAPHAM_DLIB_LUT(mvstudent),
  GRAPHAM_DLIB_LUT(nbinom), 
  GRAPHAM_DLIB_LUT(norm), 
  GRAPHAM_DLIB_LUT(pareto), 
  GRAPHAM_DLIB_LUT(poisson),
  GRAPHAM_DLIB_LUT(rayleigh),
  GRAPHAM_DLIB_LUT(student), 
  GRAPHAM_DLIB_LUT(uniform),
  GRAPHAM_DLIB_LUT(wishart), 
  {NULL, NULL, NULL}
};

/* Just a mean of the arguments */
void fmean(double *x, int Nx, 
           const double** arg, const int* len, int Nargs) { /*{{{*/
  int j, k = 0, a = 0;
  while(k<Nx && a<Nargs) {
    for (j=0; j<len[a]; j++) {
      x[k] = arg[a][j];
      k++;
      if (k>=Nx) break;
    }
    a++;
  }
  for (j=k; j<Nx; j++) {
    x[j] = 0;
  }
} /*}}}*/

void fmom2(double *x, int Nx,
            const double** arg, const int* len, int Nargs) { /*{{{*/
  int j, i, k, l, m = 0;
  for (i=0; i<Nargs; i++) {
    for (k=0; k<len[i]; k++) {
      for (j=0; j<i; j++) {
        for (l=0; l<len[j]; l++) {
          x[m] = arg[i][k]*arg[j][l];
          m++;
          if (m>=Nx) goto skip;
        }
      }
      for (l=0; l<=k; l++) {
        x[m] = arg[i][k]*arg[j][l];        
        m++;
        if (m>=Nx) goto skip;
      }
    }
  }
  skip:
  for (j=m; j<Nx; j++) x[j] = 0;
} /*}}}*/

/* Define a look up table for both C and Lua functions. */
#define GRAPHAM_FLIB_LUT(fun) {#fun, &f##fun}

const grapham_flib_struct grapham_flib[] = {
  GRAPHAM_FLIB_LUT(mean),
  GRAPHAM_FLIB_LUT(mom2),
  {NULL, NULL}
};
