/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_rand.c
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

#include <math.h>
#include <stdlib.h>
#include <sys/time.h> /* For random number generator seed */

#ifdef _HAVE_DSFMT
#include <dSFMT.h>
#endif

#include "grapham_math.h"
#include "grapham_rand.h"

/**********************************************************************
 *** Random number generation *****************************************
 **********************************************************************/

void seed_random_generator(double seed) { /*{{{*/
  unsigned long seed_;
  struct timeval tv;
  if (!isfinite(seed)) {
    gettimeofday(&tv, NULL);
    seed_ = tv.tv_sec+tv.tv_usec;
  } else {
    seed_ = (unsigned long)seed;
  }
#ifdef _HAVE_DSFMT
  init_gen_rand(seed_);
#else 
#if defined __USE_SVID || defined __USE_XOPEN
  srand48(seed_);
#else
  srand(seed_);
#endif /* defined ... */
#endif /* _HAVE_DSFMT */
} /*}}}*/

double randu(void) { /*{{{*/
#ifdef _HAVE_DSFMT
  static double array[RANDU_BUF_SZ];
  static int ind = RANDU_BUF_SZ;
  if (ind > RANDU_BUF_SZ-1) {
    fill_array_open_open(array, RANDU_BUF_SZ);
    ind = 0;
  }
  return array[ind++];
  /*return genrand_open_open();*/
#else
#if defined __USE_SVID || defined __USE_XOPEN
  return drand48();
#else
  return rand()/((double)RAND_MAX + 1);
#endif /* defined ... */
#endif /* _HAVE_DSFMT */
} /*}}}*/

double randn(void) { /*{{{*/
  /* The polar form Box-Muller */
  double u1, u2, s, tmp;
  static double z;
  static int used = 1;
  if (used) {
    do {
      u1 = 2*randu()-1;
      u2 = 2*randu()-1;
      s = u1*u1 + u2*u2;
    } while (s >= 1 || s == 0);
    tmp = sqrt(-2*log(s)/s);
    z = u1*tmp;
    used = 0;
    return u2*tmp;
  } else {
    used = 1;
    return z;
  }
} /*}}}*/

void rand_gaussian(double* x, const int d) { /*{{{*/
  int k;
  for (k=0; k<d; k++) x[k] = randn();
} /*}}}*/

void rand_laplace(double* x, const int d) { /*{{{*/
  double u;
  int k;
  for (k=0; k<d; k++) {
    u = randu();
    if (u != 0.0) {
      u = u - 0.5;
    }
    x[k] = 1.0/sqrt(2.0)*SIGNUM(u)*log(1.0-2.0*fabs(u));
  }
} /*}}}*/

void rand_uniform(double* x, const int d) { /*{{{*/
  int k;
  for (k=0; k<d; k++) {
    x[k] = sqrt(12.0)*(randu()-0.5);
  }
} /*}}}*/

void rand_student(double* x, const int d) { /*{{{*/
  int k;
  double u = randn(); u = (u != 0.0) ? u : 1.0;
  for (k=0; k<d; k++) x[k] = randn()/u;
} /*}}}*/

void mvrand(const double* L, const double sc, const double* u, double *x, const int N,
             const int full_chol) { /*{{{*/
  int k,j;
  double tmp;
  if (!full_chol) {
    /* Diagonal, so we have simply: */
    for (k=0; k<N; k++) x[k] = u[k]*sc*L[IND(k,k,N)];
  } else {
    /* In-place multiply with L^T */
    for (k=N-1; k>=0; k--) {
      tmp = 0;
      for (j=0; j<=k; j++) {
        tmp += L[IND(j,k,N)]*u[j];
      }
      x[k] = sc*tmp;
    }
  }
}

/*}}}*/

double gaussian_ratio(const double* u1, const double* u2, 
                      const double sc, const int N) { /*{{{*/
  int k;
  double p=0;
  for (k=0; k<N; k++) p -= (sc*u2[k]-u1[k])*(sc*u2[k]-u1[k]);
  for (k=0; k<N; k++) p += u1[k]*u1[k];
  return exp(.5*p);
} /*}}}*/

double uniform_ratio(const double* u1, const double* u2, 
                     const double sc, const int N) { /*{{{*/
  int j;
  for (j=0; j<N; j++) {
    if (fabs(sc*u2[j]-u1[j]) > sqrt(3.0)) {
      return 0.0;
    }
  }
  return 1.0;
} /*}}}*/

double student_ratio(const double* u1, const double* u2, 
                     const double sc, const int N) { /*{{{*/
  int k;
  double p, s;
  s = 0.0; for (k=0; k<N; k++) s += (sc*u2[k]-u1[k])*(sc*u2[k]-u1[k]);
  p = 1.0 + s;
  s = 0.0; for (k=0; k<N; k++) s += u1[k]*u1[k];
  p /= 1.0 + s;
  return pow(p, -((double)N+1.0)/2.0);
} /*}}}*/

double laplace_ratio(const double* u1, const double* u2, 
                     const double sc, const int N) { /*{{{*/
  int k;
  double p=0;
  for (k=0; k<N; k++) p -= fabs(sc*u2[k]-u1[k]);
  for (k=0; k<N; k++) p += fabs(u1[k]);
  return exp(sqrt(2)*p);
} /*}}}*/

int randint(int n) { /*{{{*/
  /* Not a particularly reliable one, especially if n is large! */
  return MIN(n,floor(randu()*(n+1)));
}
/*}}}*/

