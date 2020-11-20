/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_math.c
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

#include "grapham_math.h"

/**********************************************************************
 *** Linear algebra ***************************************************
 **********************************************************************/

int chol_linpack(double *A, int N) { /*{{{*/
  double* tmp = (double *)malloc(sizeof(double)*N);
  int info, zero=0;
  dchdc_(A, &N, &N, tmp, NULL, &zero, &info);
  free(tmp);
  if (info<N) {
    return -1;
  } else {
    return 0;
  }
}
/*}}}*/

void chol_update_linpack(double* L, double alpha, 
                        double* x, double* tmp, int N) { /*{{{*/
  int zero=0, k; 
  double *s = tmp, *c = &(tmp[N]);
  double salpha = sqrt(alpha);
  for (k=0; k<N; k++) {
    x[k] *= salpha;
  }
  dchud_(L, &N, &N, x, NULL, &N, &zero, NULL, NULL, s, c);
}
/*}}}*/

int chol_downdate_linpack(double* L, double alpha, 
                        double* x, double* tmp, int N) { /*{{{*/
  int info=0;
  int zero=0, k; 
  double *s = tmp, *c = &(tmp[N]);
  double salpha = sqrt(alpha);
  for (k=0; k<N; k++) {
    x[k] *= salpha;
  }
  dchdd_(L, &N, &N, x, NULL, &N, &zero, NULL, NULL, s, c, &info);
  return(info);
}
/*}}}*/

void tril_solve_inplace(const double* c, double* x, const int N) { /*{{{*/
  double rho;
  int i,j;
  for (i=0; i<N; i++) {
    rho = 0;
    for (j=0; j<=i-1; j++) {
      rho += c[IND(j,i,N)]*x[j];
    }
    x[i] = (x[i]-rho)/c[IND(i,i,N)];
  }
}

/*}}}*/

/* Compute the inverse (and return the log-determinant) of an 
 * uppper-triangular matrix. */
double triu_inv_linpack(double *v, int p) { /*{{{*/
  double det[2];
  int job = 111;
  int info;
  dtrdi_(v, &p, &p, det, &job, &info);
  if (info != 0) {
    return GRAPHAM_NINF;
  } else {
    return LN10*det[1]+log(det[0]);
  }
} /*}}}*/

/* Compute the log-determinant of a triangular matrix. */
double tri_det_linpack(double *v, int p) { /*{{{*/
  double det[2];
  int job = 100;
  dtrdi_(v, &p, &p, det, &job, NULL);
  return LN10*det[1]+log(det[0]);
} /*}}}*/

void scale_triu(double* m, const double s, const int N) { /*{{{*/
  int i,j;
  for (i=0; i<N; i++) {
    for (j=0; j<=i; j++) {
      m[IND(j,i,N)] *= s;
    }
  }
}
/*}}}*/

void sub_vector(double* z, const double *x, const double* y, const int N) { /*{{{*/
  int k;
  for (k=0; k<N; k++) z[k] = x[k] - y[k];
}
/*}}}*/

void mac_vector(double *z, const double *x, const double s, const double* y, int N) { /*{{{*/
  int k;
  for (k=0; k<N; k++) z[k] = x[k] + s*y[k];
}
/*}}}*/

void set_identity_matrix(double *m, const int N) { /*{{{*/
  int i,j;
  for (j=0; j<N; j++) {
    for (i=0; i<N; i++) {
      m[N*j+i] = (j == i)?1.0:0.0;
    }
  }
}

/*}}}*/

/**********************************************************************
 *** Continuous multidimensional distributions ************************
 **********************************************************************/

double d_mvnorm(const double* x, const double* m, double* v, 
                    int N) { /*{{{*/
  double* y, p;
  int tmp;

  y = (double *)malloc(sizeof(double)*N);
  
  tmp = chol_linpack(v, N);
  if (tmp < 0) return(GRAPHAM_NINF);
  
  p = d_mvnorm_chol_noalloc(x, m, v, y, N);
  free(y);
  return(p);
}

/*}}}*/

double d_mvnorm_chol_noalloc(const double* x, const double* m, 
                            const double* c, double* y, int N) { /*{{{*/
  int k;
  double p = -.5*LN2PI*N;
  for (k=0; k<N; k++) {
    y[k] = x[k]-m[k];
  }
  tril_solve_inplace(c, y, N);
  for (k=0; k<N; k++) {
    p -= .5*y[k]*y[k] + log(c[IND(k,k,N)]);
  }
  return(p);
}

/*}}}*/

double d_mvnorm_chol(const double* x, const double* m, const double* c, 
                    int N) { /*{{{*/
  double* y;
  double p;
  y = (double *)malloc(sizeof(double)*N);
  p = d_mvnorm_chol_noalloc(x, m, c, y, N);
  free(y);
  return(p);
}

/*}}}*/

double d_mvnorm_diag(const double* x, const double* m, const double* v, 
                    int N) { /*{{{*/
  int k;
  double p = 0;
  
  for (k=0; k<N; k++) {
    if (v[k] <= 0) return(GRAPHAM_NINF);
    p += d_norm(x[k], m[k], v[k]);
  }
  return(p);
}

/*}}}*/

double lnmvgamma(double x, double N) { /*{{{*/
  int j;
  double y, r = .25*N*(N-1.0)*LNPI;
  
  for (j=0; j<N; j++) {
    y = x+(1.0-(double)j)/2.0;
    r += dlgama_(&y);
  }
  return r;
} /*}}}*/

double d_mvstudent(double* x, double* mu, double* v, double n, int p) { /*{{{*/
  int k;
  double vdet = 0, d = 0;
  double n2 = .5*n, np2 = .5*(n+p);
  int cond;
  
  if (n<=0 || p<=0) return GRAPHAM_NINF;
  
  cond = chol_linpack(v, p);
  if (cond < 0) {
    return GRAPHAM_NINF;
  } else {
    for (k=0; k<p; k++) {
      x[k] -= mu[k];
    }
    tril_solve_inplace(v, x, p);
    for (k=0; k<p; k++) {
      vdet += log(v[IND(k,k,p)]);
      d += x[k]*x[k];
    }
    return dlgama_(&np2) - dlgama_(&n2) 
    - .5*(p*log(n) + p*log(PI) + vdet + (n+p)*log(1.0+1.0/n*d));
  }
}
/*}}}*/

double d_wishart(double* x, double* v, const double n, int p) { /*{{{*/
  double lndetv, lndetx, tr, vkj, xjk;
  int k, j, l, condx, condv;
  
  if (n<=0 || p<=0) return GRAPHAM_NINF;
  
  condx = chol_linpack(x, p);
  condv = chol_linpack(v, p);
  
  if (condx < 0 || condv < 0) {
    return GRAPHAM_NINF;
  } else {
    lndetv = triu_inv_linpack(v, p);
    lndetx = tri_det_linpack(x, p);
    if (!isfinite(lndetv) || !isfinite(lndetx)) return GRAPHAM_NINF;
    tr = 0;
    for (k=0; k<p; k++) {
      for (j=0; j<p; j++) {
        vkj = 0; xjk = 0;
        for (l=0; l<=k && l<=j; l++) {
          vkj += v[IND(l,k,p)]*v[IND(l,j,p)];
          xjk += x[IND(l,k,p)]*x[IND(l,j,p)];
        }
        tr += vkj*xjk;
      }
    }
    return lndetx*(n-p-1.0) - lndetv*n - lnmvgamma(.5*n, p)
    + .5*(-n*p*LN2-tr);
  }
} /*}}}*/

/**********************************************************************
 *** Continuous one-dimensional distributions *************************
 **********************************************************************/

double log_beta(double a, double b) { /*{{{*/
  double a_b = a+b;
  return dlgama_(&a) + dlgama_(&b) - dlgama_(&a_b);
}
/*}}}*/

double d_beta(double x, double alpha, double beta) { /*{{{*/
  if (alpha <= 0 || beta <= 0 || x <= 0 || x >= 1) {
    return GRAPHAM_NINF;
  } else {
    return (alpha-1)*log(x) + (beta-1)*log(1-x) 
           -log_beta(alpha, beta);
  }
}

/*}}}*/

double d_chi2(double x, double k) { /*{{{*/
  double k_half;
  if (x<0 || k<= 0) {
    return GRAPHAM_NINF;
  } else {
    k_half = .5*k;
    return k_half*log(0.5) - dlgama_(&k_half) + (k_half-1)*log(x) - .5*x;
  }
} /*}}}*/

double d_cauchy(double x, double x0, double gamma) { /*{{{*/
  double xdiff = x-x0, p = -1.14472988584940; /* -log(pi) */
  if (gamma<=0) {
    p = GRAPHAM_NINF;
  } else {
    /*double xnorm = (x-x0)/gamma;
     return -log(PI*gamma*(1+xnorm*xnorm));*/
    p += log(gamma) - log(gamma*gamma + xdiff*xdiff);
  }
  return p;
} /*}}}*/

double d_exp(double x, double s) { /*{{{*/
  if (x<0 || s<=0) {
    return(GRAPHAM_NINF);
  } else {
    return(-x/s-log(s));
  }
}

/*}}}*/

double d_erlang(double x, double k, double lambda) { /*{{{*/
  double k1;
  if (k<=0 || lambda<=0 || x<0) {
    return GRAPHAM_NINF;
  } else {
    k1 = k+1;
    return k*log(lambda) + (k-1)*log(x) - lambda*x - dlgama_(&k1);
  }
}
/*}}}*/

double d_fisher(double x, double d1, double d2) { /*{{{*/
  if (x<0 || d1 <= 0 || d2 <= 0) {
    return GRAPHAM_NINF;
  } else {
    return .5*( d1*log(d1*x) + d2*log(d2) - (d1+d2)*log(d1*x+d2) )
           -log(x) -log_beta(d1/2, d2/2);
  }
} /*}}}*/

double d_gamma(double x, double k, double theta) { /*{{{*/
  if (x<0 || theta <= 0 || k <= 0) {
    return GRAPHAM_NINF;
  } else {
    return (k-1)*log(x) - x/theta -k*log(theta) - dlgama_(&k);
  }
}

/*}}}*/

double d_gumbel(double x, double mu, double beta) { /*{{{*/
  double z;
  if (beta<=0) {
    return GRAPHAM_NINF;
  } else {
    z = -(x-mu)/beta;
    return z-exp(z)-log(beta);
  }
} /*}}}*/

double d_invchi2(double x, double nu) { /*{{{*/
  return d_invgamma(x, nu/2, 1/2);
}
/*}}}*/

double d_invgamma(double x, double alpha, double beta) { /*{{{*/
  if (x<=0 || alpha <= 0 || beta <= 0) {
    return GRAPHAM_NINF;
  } else {
    return alpha*log(beta) - dlgama_(&alpha) - (alpha+1)*log(x)
           - beta/x;
  }
}
/*}}}*/

double d_laplace(double x, double mu, double b) { /*{{{*/
  if (b<=0) {
    return GRAPHAM_NINF;
  } else {
    return -abs(x-mu)/b-log(2*b);
  }
}

/*}}}*/

double d_levy(double x, double c) { /*{{{*/
  if (x <= 0 || c <= 0) {
    return GRAPHAM_NINF;
  } else {
    return .5*(log(c/(2*PI)) -.5*c/x - 3*log(x));
  }
} /*}}}*/

double d_logistic(double x, double mu, double s) { /*{{{*/
  double u;
  if (s<=0) {
    return GRAPHAM_NINF;
  } else {
    u = -(x-mu)/s;
    return u - log(s) - 2*log(1+exp(u));
  }
} /*}}}*/

double d_lognorm(double x, double m, double v) { /*{{{*/
  double d, lnx;
  if (x<=0 || v<=0) {
    return GRAPHAM_NINF;
  } else {
    lnx = log(x);
    d = lnx-m;
    return -lnx -.5*(LN2PI + d*d/v + log(v));
  }
} /*}}}*/

double d_norm(double x, double m, double v) { /*{{{*/
  double d;
  if (v<=0) return(GRAPHAM_NINF);
  d = x-m;
  return -.5*(LN2PI + d*d/v + log(v));
}

/*}}}*/

double d_pareto(double x, double xm, double k) { /*{{{*/
  if (xm<=0 || x<xm || k<=0) {
    return GRAPHAM_NINF;
  } else {
    return -(k+1)*log(x) + log(k) + k*log(xm);
  }
}
/*}}}*/

double d_rayleigh(double x, double s) { /*{{{*/
  double s2;
  if (s<=0) {
    return GRAPHAM_NINF;
  } else {
    s2 = s*s;
    return log(x) - x*x/(2*s2) - 2*log(s2);
  }
} /*}}}*/

double d_student(double x, double nu) { /*{{{*/
  double nu_1p2, nu_2;
  if (nu<=0) {
    return GRAPHAM_NINF;
  } else {
    nu_1p2 = (nu+1)/2; nu_2 = nu/2;
    return -0.57236494292470008707 /* log(1/sqrt(pi)) */
           + dlgama_(&nu_1p2)-dlgama_(&nu_2) - .5*log(nu)
           -(nu+1)/2 * log(1 + x*x/nu);
  }
} /*}}}*/

/**********************************************************************
 *** Discrete one-dimensional distributions ***************************
 **********************************************************************/

double log_nchoosek(double n, double k) { /*{{{*/
  double n_ = n+1, k_ = k+1, nk_ = n-k+1;
  if (k==0) {
    return 0;
  } else {
    return dlgama_(&n_) - dlgama_(&k_) - dlgama_(&nk_);
  }
} /*}}}*/

double d_binom(double k, double n, double p) { /*{{{*/
  if (p<0 || p>1 || n<0 || k<0 || k>n) return GRAPHAM_NINF;
  return log_nchoosek(n, k) + k*log(p) + (n-k)*log(1-p);
}

/*}}}*/

double d_nbinom(double k, double r, double p) { /*{{{*/
  double r_k = r + k;
  double k_ = k+1;
  if (k<0 || r<=0 || p<=0 || p>=1) {
    return GRAPHAM_NINF;
  } else {
    return dlgama_(&r_k) - dlgama_(&k_) - dlgama_(&r)
           + r*log(p) + k*log(1-p);
  }
}

/*}}}*/

double d_poisson(double k, double lambda) { /*{{{*/
  double k_ = k+1;
  if (k < 0 || lambda <= 0) {
    return GRAPHAM_NINF;
  } else {
    return -dlgama_(&k_) - lambda + k*log(lambda);
  }
}

/*}}}*/
