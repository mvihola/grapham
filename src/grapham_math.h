/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_math.h
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

/* Some definitions of constants (if not defined in math.h). */

#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif
#define LN10 2.302585092994045901093613792
#define LN2 6.9314718055994528622676398299518041312695e-01
#define LNPI 1.1447298858494001638774761886452324688435e+00
#define LN2PI 1.8378770664093453390819377091247588396072e+00

#ifndef INFINITY
/* We assume IEEE floating point arithmetic, in which case 
 * these yield infinity and negative infinity. */
#define GRAPHAM_INF (+1.0/+0.0)
#define GRAPHAM_NINF (-1.0/+0.0)
#else
#define GRAPHAM_INF INFINITY
#define GRAPHAM_NINF -INFINITY
#endif

#define MAX(a,b) ((a>b) ? a : b)
#define MIN(a,b) ((a<b) ? a : b)
#define SIGNUM(a) ((a<0.0) ? -1.0 : 1.0)

/* IND(i,j,N) is the (i,j):th element of a matrix with
 * first dimension N. Note: this corresponds to the 
 * column-major order used in Fortran! */
#define IND(i,j,N) (i+N*(j))

/**********************************************************************
 *** Netlib externals *************************************************
 **********************************************************************/

/* Cholesky factorisation in netlib linpack. */
extern void dchdc_( double* a, int* lda, int* p, double* work, 
                    int* jpvt, int* job, int* info); 

/* Cholesky update in netlib linpack. */
extern void dchud_(double* r,int* ldr,int* p,double* x,double* z,
                   int* ldz,int* nz,double* y,double *rho,double* c,double *s);

/* Cholesky downdate in netlib linpack. */
extern void dchdd_(double* r,int* ldr,int* p,double* x,double* z,
                   int* ldz,int* nz,double* y,double *rho,double* c,double *s, int* info);

/* Log-gamma function of netlib specfun. */
extern double dlgama_(double* x);

/* Determinant and inverse of a triangular matrix. */
extern void dtrdi_( double *t, int* ldt, int* n, double* det, int* job, 
                    int* info);

/**********************************************************************
 *** Linear algebra ***************************************************
 **********************************************************************/

/* In-place Cholesky factorisation.
 * 
 * In:
 *   A -- Symmetric and positive definite matrix. Note: only the 
 *        upper-triangular part of A is used.
 *   d -- Dimension of A.
 *
 * Out:
 *   A -- Upper-triangular Cholesky factor of input matrix A.
 *   status -- Zero if succesful, -1 if not.
 */
int chol_linpack(double* A, int d);

/* Cholesky rank-1 update.
 * 
 * In:
 *   L     -- The Cholesky factor of A, i.e. A = L^T L.
 *   alpha -- A non-negative multiplier of the update vector.
 *   x     -- The update vector.
 *   tmp   -- Temporary storage of size (at least) twice the dimension.
 *   d     -- The dimension.
 * 
 * Out:
 *   L     -- The Cholesky factor of  A + alpha x x^T.
 */
void chol_update_linpack(double* L, double alpha, double* x,
                         double* tmp, int d);

/* Cholesky rank-1 downdate.
 * (Same syntax as with chol_update_linpack, but L on output is 
 * the Cholesky factor of  A - alpha x x^T.) */
int chol_downdate_linpack(double* L, double alpha, double* x,
                           double* tmp, int d);

/* Solve y from the linear system C^T y = x, where C is upper triangular.
 * 
 * In:
 *   C -- The upper-triangular matrix.
 *   x -- The constant vector.
 *   d -- The dimension.
 * 
 * Out:
 *   x -- The solution y.
 */
void tril_solve_inplace(const double* C, double* x, const int d);

/* Scale an upper-triangular matrix.
 * 
 * In:
 *   M -- The upper-triangular matrix.
 *   s -- the scaling factor.
 *   d -- The dimension.
 * 
 * Out:
 *   M -- The input matrix M scaled by s. 
 */
void scale_triu(double* M, const double s, const int d);

/* Subtract a vector from another.
 *
 * In:
 *   x, y -- Two vectors.
 *   d    -- The dimension.
 * 
 * Out:
 *   z -- The vector x-y.
 * 
 * The output array "z" can be either of the inputs. 
 */
void sub_vector(double *z, const double *x, const double* y, const int d);

/* Multiply and accumulate a vector.
 * 
 * In:
 *   x, y -- Two vectors.
 *   s    -- The scaling factor.
 *   d    -- The dimension.
 * 
 * Out:
 *   z    -- The vector x + s*y.
 * 
 * The output array "z" can be either of the inputs. 
 */
void mac_vector(double* z, const double *x, const double s, const double* y, 
                const int d);

/* Set to an identity matrix. 
 * 
 * In:
 *   d -- The dimension.
 * 
 * Out:
 *   M -- The d-by-d identity matrix.
 */
void set_identity_matrix(double *M, const int d);

/**********************************************************************
 *** Continuous multidimensional distributions ************************
 **********************************************************************/

/* Multivariate normal with mean vector m and positive definite covriance
 * matrix v. Versions: */
/* 1) Diagonal covariance: 
 *    v is a vector containing the positive diagonal elements. */
double d_mvnorm_diag(const double* x, const double* m, const double* v, 
                    int N);
/* 2) Cholesky factor of the covariance: 
 *    v is upper-triangular Cholesky factor of the covariance matrix. */
double d_mvnorm_chol(const double* x, const double* m, const double* L, 
                    int N);
double d_mvnorm_chol_noalloc(const double* x, const double* m, const double* L, 
                           double* tmp, int N);
/* 3) Full covariance:
 *    The upper-triangular elements of v contain the covariance matrix. */
double d_mvnorm(const double* x, const double* m, double* v, 
                    int N);

/* Multivariate student's t-distribution with parameter
 * vector mu, positive definite v, and real p>0. */
double d_mvstudent(double* x, double* mu, double* v, double p, int N);

/* Wishart distribution */
double d_wishart(double* x, double* v, const double n, int p);

/**********************************************************************
 *** Continuous one-dimensional distributions *************************
 **********************************************************************/

/* Beta with shapes alpha>0 and beta>0. */
double d_beta(double x, double alpha, double beta);

/* Chi-squared with k>0 degrees of freedom. */
double d_chi2(double x, double k);

/* Cauchy-Lorentz with location x0 and scale gamma>0. */
double d_cauchy(double x, double x0, double gamma);

/* Erlang with shape k>0 and rate lambda>0. */
double d_erlang(double x, double k, double lambda);

/* Exponential with scale (inverse rate) b>0. */
double d_exp(double x, double b);

/* F-distribution ("Fisher-Snedecor") with d1>0 and d2>0 degrees
 * of freedom. */
double d_fisher(double x, double d1, double d2);

/* Gamma with shape k>0 and scale theta>0. */
double d_gamma(double x, double k, double theta);

/* Gumbel (Fisher-Tippett), with location mu and scale beta>0. */
double d_gumbel(double x, double mu, double beta);

/* Inverse-gamma with shape alpha>0 and scale beta>0. */
double d_invgamma(double x, double alpha, double beta);

/* Inverse-chi-squared with nu>0 degrees of freedom. */
double d_invchi2(double x, double nu);

/* Laplace with mean mu and scale (inverse rate) b>0. */
double d_laplace(double x, double mu, double b);

/* Lévy with scale c>0. */
double d_levy(double x, double c);

/* Logistic with location mu and scale s>0. */
double d_logistic(double x, double mu, double s);

/* Log-normal with mean m and variance v>0 of log(x). */
double d_lognorm(double x, double mu, double v);

/* Normal with mean m and variance v>0. */
double d_norm(double x, double m, double v);

/* Pareto with location xm>0 and shape k>0. */
double d_pareto(double x, double xm, double k);

/* Rayleigh with scale s>0. */
double d_rayleigh(double x, double s);

/* Student's t with nu>0 degrees of freedom. */
double d_student(double x, double nu);

/**********************************************************************
 *** Discrete one-dimensional distributions ***************************
 **********************************************************************/

/* Binomial with n>=0 trials and success probability 0<=p<=1. */
double d_binom(double k, double n, double p);

/* Negative binomial with parameters r>0 and 0<p<1. */
double d_nbinom(double k, double r, double p);

/* Poisson with rate lambda>0. */
double d_poisson(double x, double lambda);
