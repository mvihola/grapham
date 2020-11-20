/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_rand.h
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

/**********************************************************************
 *** Random number generation *****************************************
 **********************************************************************/

/* The buffer for increasing dSFMT calling efficiency. */
#define RANDU_BUF_SZ 65536

/* Seed the random number generator.
 * 
 * In:
 *   seed -- The seed of the random number generator.
 */
void seed_random_generator(double seed);

/* Uniform [0,1) random variables with drand48() or with dSFMT. */
double randu(void);

/* Standard Gaussian N(0,1) random variables with Box-Muller transform and
 * randu(). */
double randn(void);

/* Multivariate N(0,I_d) Gaussian random variates. 
 *
 * In:
 *   x -- The array where the random vector is written.
 *   d -- The dimension of the vector.
 */
void rand_gaussian(double* x, const int d);

/* Vector of independent uniform (-sqrt(3),sqrt(3)) random variates.
 * (See rand_gaussian for details.) */
void rand_uniform(double* x, const int d);

/* Multivariate Student with one degree of freedom (a.k.a. multivariate 
 * Cauchy random variates. (See rand_gaussian for details.) */
void rand_student(double* x, const int d);

/* Vector of independent Laplace random variates with mean 0, variance 1.
 * (See rand_gaussian for details.) */
void rand_laplace(double* x, const int d);

/* Compute a scaled random vector x = sc L^T u, where sc is a scalar 
 * and L is an upper-triangular matrix.  
 *
 * In:
 *   L  -- The upper-triangular matrix.
 *   sc -- The scaling factor.
 *   u  -- A vector of dimension d.
 *   x  -- The array where the random vector is written.
 *   d  -- The dimension of the vector.
 *   full_chol -- If zero, the matrix L is assumed diagonal. 
 */
void mvrand(const double* L, const double sc, const double* u, double* x, const int d, 
             const int full_chol);

/* The ratio of the first and second stage proposal in DRAM,
 * when using a Gaussian proposal.
 * 
 * In:
 *   u1 -- First-stage proposal vector *before* scaling.
 *   u2 -- Second-stage proposal vector *before* scaling.
 *   sc -- The DR scaling factor; second stage proposal is the first
 *         stage proposal scaled by sc.
 *   d  -- The dimension.
 * 
 * Out: q(y1-y2)/q(y2-x), where  y1 = r L^T u1 + x  and  
 *      y2 = sc r L^T u2 + x  with some scalar r>0 and 
 *      a non-singular matrix L^T.
 */
double gaussian_ratio(const double* u1, const double* u2, 
                      const double sc, const int d);

/* The ratio of the first and second stage proposal in DRAM,
 * when using uniform proposal. (See gaussian_ratio for details) */
double uniform_ratio(const double* u1, const double* u2, 
                     const double sc, const int d);

/* The ratio of the first and second stage proposal in DRAM,
 * when using student proposal. (See gaussian_ratio for details) */
double student_ratio(const double* u1, const double* u2, 
                     const double sc, const int d);

/* The ratio of the first and second stage proposal in DRAM,
 * when using laplace proposal. (See gaussian_ratio for details) */
double laplace_ratio(const double* u1, const double* u2, 
                     const double sc, const int d);

/* Random integers from the interval [0,n] */
int randint(int n);
