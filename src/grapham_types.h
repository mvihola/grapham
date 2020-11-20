/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_types.h
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

#ifndef _GRAPHAM_TYPES_H
#define _GRAPHAM_TYPES_H

#ifndef _GRAPHAM_VERSION
#define _GRAPHAM_VERSION ???
#endif

/* Constants */ /*{{{*/
#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

/* Default scaling of a block */
#define DEFAULT_SCALING(d) ((2.38)/sqrt(d))

#define USAGE_TEXT \
"This program comes with ABSOLUTELY NO WARRANTY. This is free\n" \
"software, and you are welcome to redistribute it under certain\n" \
"conditions; see the included file COPYING for details.\n\n" \
"Usage: %s [-v|-vv|-q] [-e \"lua code\"] [model_file(s)]\n"

#define QUOTE(s) QQUOTE(s)
#define QQUOTE(s) #s
#define VERSION_TEXT \
"Grapham " QUOTE(_GRAPHAM_VERSION) "  Copyright (C) 2009-2013  Matti Vihola.\n" 

#define LUA_INIT_SCRIPT "init/grapham_init.lua"
#define LUA_PREPROCESS_SCRIPT "init/grapham_preprocess.lua"

#define OUTFILE_NUMFORMAT "%.16e"
#define OUTFILE_FIELD_SEP ","
#define OUTFILE_RECORD_SEP "\n"
#define OUTCFG_NUMFORMAT "%.16e,"

#define OUTFMT_ASCII  "ascii"
#define OUTFMT_BINARY  "bin"

#define ALG_AM_NAME "am"
#define ALG_ASCM_NAME "asm"
#define ALG_AMS_NAME "aswam"
#define ALG_RBAM_NAME "rbam"
#define ALG_RBAMS_NAME "rbaswam"
#define ALG_ASHM_NAME "ram"
#define ALG_METROPOLIS_NAME "metropolis"

#define BLOCKING_SC_NAME "sc"
#define BLOCKING_NODE_NAME "node"
#define BLOCKING_FULL_NAME "full"

#define TYPE_NUMBER_NAME "number"
#define TYPE_VECTOR_NAME "vector"
#define TYPE_CUSTOM_NAME "custom"
#define TYPE_MATRIX_NAME "matrix"

#define KIND_REAL_NAME "real"
#define KIND_INTEGER_NAME "integer"
#define KIND_SYMMETRIC_NAME "symmetric"

#define INIT_GREEDY_NAME "greedy"
#define INIT_FREEZE_NAME "freeze"
#define INIT_TRAD_NAME "trad"

#define CONST_VNAME "const"

#define MODEL_VNAME  "model"
#define MODEL_TYPE_VNAME  "type"
#define MODEL_KIND_VNAME "kind"
#define MODEL_DIM_VNAME  "dim"
#define MODEL_DENSITY_VNAME  "density"
#define MODEL_PARENTS_VNAME  "parents"
#define MODEL_INIT_VNAME  "init_val"
#define MODEL_INSTANTIATED_VNAME "instantiated"
#define MODEL_LIMITS_VNAME "limits"

#define INST_VNAME  "data"

#define PROPOSAL_GAUSSIAN "norm"
#define PROPOSAL_LAPLACE "laplace"
#define PROPOSAL_UNIFORM "unif"
#define PROPOSAL_STUDENT "student"

#define PARA_VNAME  "para"
#define PARA_NITER_VNAME  "niter"
#define PARA_NBURN_VNAME  "nburn"
#define PARA_NTHIN_VNAME  "nthin"
#define PARA_ACC1_VNAME  "acc_opt"
#define PARA_ACC2_VNAME  "acc_opt2"
#define PARA_SEED_VNAME "seed"
#define PARA_ALG_VNAME  "algorithm"
#define PARA_BLOCKING_VNAME "blocking"
#define PARA_BLOCKS_VNAME  "blocks"
#define PARA_BLOCKS_CHOL_VNAME  "blocks_chol"
#define PARA_BLOCKS_CHOL_SC_VNAME  "blocks_chol_sc"
#define PARA_BLOCKS_SC_VNAME  "blocks_scaling"
#define PARA_OUTF_VNAME  "outfile"
#define PARA_OUTFMT_VNAME  "outfmt"
#define PARA_OUTVARS_VNAME  "outvars"
#define PARA_CLOSE_VNAME  "close_hook"
#define PARA_CUSTOM_SCALING_FUN_VNAME "scaling_adapt"
#define PARA_OUTADAPT_VNAME "adapt_outfile"
#define PARA_OUTADAPTFMT_VNAME "adapt_outfmt"
#define PARA_OUTCFG_VNAME "outcfg"
#define PARA_DR_SCALING_VNAME "dr"
#define PARA_CLIB_VNAME "clib"
#define PARA_ADAPT_WEIGHT_VNAME "adapt_weight"
#define PARA_ADAPT_WEIGHT_SC_VNAME "adapt_weight_sc"
#define PARA_INIT_VNAME "init"
#define PARA_MIX_WEIGHT_VNAME "p_mix"
#define PARA_PROPOSAL_VNAME "proposal"
#define PARA_RANDOM_SCAN_VNAME "random_scan"

#define FUNC_VNAME  "functional"
#define FUNC_DIM_VNAME "dim"
#define FUNC_NAME_VNAME "name"
#define FUNC_ARGS_VNAME "args"

#define DEFAULT_NITER    ((int)10e3)
#define DEFAULT_NBURN    (0)
#define DEFAULT_NTHIN    (1)
#define DEFAULT_ACC_OPT1 (0.44)
#define DEFAULT_ACC_OPT2 (0.234)
#define DEFAULT_ALG      (ALG_AM)
#define DEFAULT_BLOCKING (BLOCKING_NODE)
#define DEFAULT_INIT     (INIT_GREEDY)
#define DEFAULT_DR_SCALING (0.0)
#define DEFAULT_KIND     (REAL_KIND)
#define DEFAULT_PROPOSAL (PROPOSAL_GAUSSIAN)
#define DEFAULT_RANDOM_SCAN (0)
#define DEFAULT_ADAPT_WEIGHT (1.0)
#define DEFAULT_ADAPT_WEIGHT_SC (0.66)

/*}}}*/

/* Types */ /*{{{*/
typedef enum { FALSE = 0, TRUE } boolean;

typedef enum { /*{{{*/
  ALG_AM,       /* Adaptive Metropolis */
  ALG_AMS,      /* Adaptive Metropolis with adaptive scaling */
  ALG_RBAM,     /* Rao-Blackwellised Adaptive Metropolis */
  ALG_RBAMS,    /* RBAM with adaptive scaling */
  ALG_ASCM,     /* Adaptive scaling Metropolis */
  ALG_ASHM,     /* Adaptive shape Metropolis */
  ALG_METROPOLIS /* Metropolis (one sampling block) */
  /*}}}*/
} algorithm;

typedef enum { /*{{{*/
  INIT_GREEDY,
  INIT_FREEZE,
  INIT_TRAD
  /*}}}*/
} init;

typedef enum { /*{{{*/
  BLOCKING_SC,
  BLOCKING_NODE,
  BLOCKING_FULL
/*}}}*/
} blocking;

typedef enum { /*{{{*/
  UNKNOWN_TYPE, /* Not yet determined */
  NUMBER_TYPE,  /* a Lua number */
  VECTOR_TYPE,  /* a Lua table */
#ifdef _HAVE_NUMLUA
  MATRIX_TYPE,  /* Numlua matrix */
#endif
  CUSTOM_TYPE   /* Anything indexed as a Lua table */
/*}}}*/
} Var_type;

typedef enum { /*{{{*/
  REAL_KIND,
  SYMMETRIC_KIND,
  INTEGER_KIND
/*}}}*/
} Var_kind;

typedef struct Node_el { /*{{{*/
  Var_type type;     /* The type of the variable */
  Var_kind kind;     /* The kind of the variable */
  char *vname;       /* The name of the variable */
  int dim;           /* The total dimension of the variable */
  int Ndims;         /* Number of dimensions of the variable */
  int* size;         /* Size of each dimension */
  double *value;     /* The value of the variable */
  int lua_value;     /* The Lua registry index of the value */
  int density;       /* The Lua registry index of the log-density function */
  double (*cdensity)(double **, int *, int); /* The C-density, if exists */
  boolean is_cdensity; /* Whether there is a C-density. */
  double **cdensity_args;
  int *cdensity_lens;
  double density_value; /* The current value of the node density */
  double new_density_value; /* The density value given the proposal */
  int Nparents;      /* Number of parents */
  struct Node_el** parents; /* List of pointers to the parent nodes */
  int Nchildren;     /* Number of children */
  struct Node_el** children; /* List of pointers to the child nodes */
  boolean instantiated; /* True, if the variable is instantiated */
  boolean* in_block; /* For each element, information whether the element is in
                      * a sampling block */
  boolean connected; /* This is used by the graph acyclicity & connectivity
                      * checking algorithms */
  boolean censored;  /* Whether there is censoring */
  double limits[2];  /* The censoring limits */
/*}}}*/
} Node;

typedef struct { /*{{{*/
  Node* node; /* The node having the variable information */
  int index;  /* The index, i.e. we use node->value[index] */
  /*}}}*/
} Component;

typedef struct { /*{{{*/
  int dim;           /* The dimension of the block */
  Component* components; /* Variable information for each component */
  int Nvars;         /* Number of variables in the block */
  Node** variables;  /* List of pointers to the variable nodes in the block */
  int Nchildren;     /* Number of children of the whole block. */
  Node** children;   /* Pointers to the child nodes. */
  double *old_value; /* Temporary storage (of length dim) of the old and the */
  double *proposal_value; /* Temporary storage of proposed value during the proposal step. */
  double *rand_value; /* Temporary unscaled proposal increment value. */
  double *value;     /* Pointer to the "current" value in the proposal step
                        (either old_value or proposal_value). */
  double *tmp_value;  /* Temporary values; Length: dim*3 */
  double scaling;    /* The adapted scaling parameter. */
  double init_scaling;    /* The initial (fixed) scaling parameter. */
  double *adapt_mean; /* The adapted mean parameter. */
  double *init_chol;      /* The initial (fixed) Cholesky factor parameter. */
  double *adapt_chol; /* The adapted Cholesky factor parameter. */
  unsigned int accepted; /* How many times a proposal of this block is accepted */
  unsigned int rejected; /* How many times a proposal of this block is rejected */
  /*}}}*/
} Block;

typedef struct { /*{{{*/
  int N;         /* The number of nodes in the model. */
  Node* nodes;   /* An array of N nodes in the model. */
  int dim;       /* The total dimension of the model. */
  Block** blocks; /* The sampling blocks */
  int Nblocks;   /* Total number of sampling blocks */
  /*}}}*/
} Model;

typedef struct Output_ { /*{{{*/
    double** vars;   /* The saved components */
    int Nvars;       /* Number of saved components */
    FILE* file;      /* The output file handle */
    boolean write_file; /* A flag whether output is requested */
    void (*add_record_outfile)(struct Output_*); /* The function writing one
                                                 * output record */
/*}}}*/
} Output;

typedef struct Parameters_ { /*{{{*/
  unsigned long niter; /* Number of actual iterations */
  unsigned long nburn; /* Number of burn-in iterations */
  unsigned long nthin; /* Number of samples in thinning */
  double acc_opt1;     /* "Optimal" acceptance rate in 1D-case */
  double acc_opt2;     /* "Optimal" acceptance rate in >=2D-case */
  double dr_scaling;   /* The scaling of the second stage of DR. */
  double seed;         /* The seed supplied to the random number generator */
  algorithm alg;       /* The selected sampling algorithm */
  blocking blocking;   /* The selected blocking strategy */
  init init;           /* The selected initialisation strategy */
  int close_fun;       /* A function (Lua ref) to be called when finished. */
  boolean is_close_fun;/* Whether there is a close_fun */
  int custom_scaling_fun;
  boolean is_custom_scaling_fun;
  int adapt_weight_fun; /* AM related weight */
  boolean is_adapt_weight_fun;
  double adapt_weight_exp;
  int adapt_weight_sc_fun; /* Scaling weight */
  boolean is_adapt_weight_sc_fun;
  double adapt_weight_sc_exp;
  int mix_weight_fun;
  boolean is_mix_weight_fun;
  boolean random_scan;
  void *clib;          /* The handle of the C-library (if exists) */
  FILE* outcfg;        /* Output configuration file handle (or NULL) */
  /* A function adapting the scaling parameter */
  double (*scaling_adapt)(lua_State *, double, double, double, int, int, 
                          struct Parameters_ *);
  /* A function updating the (mean, Cholesky factor) pair of a block. */
  void (*mean_chol_update)(Block*, struct Parameters_ *, double, double, double);
  /* A function for one accept-reject step of a block */
  double (*accept_reject)(lua_State*, Block*, struct Parameters_ *, double*, double);
  /* Function generating an unscaled proposal vector. */
  void   (*rand)(double*, const int);
  /* The proposal ratio of the first and second stage proposals. */
  double (*dram_proposal_ratio)(const double *, const double*, const double, const int);
  Output output;
  Output adapt_output;
  /*}}}*/
} Parameters;

typedef struct { /*{{{*/
  double *mean;
  int mean_N;
  int fun;
  void (*cfun)(double*, int, const double **, const int *, int);
  boolean is_cfun;
  const double** cfun_args;
  double* cfun_tmp;
  int* cfun_lengths;
  int cfun_Nargs;
  Var_type type;
  boolean evaluate;
  unsigned long called;
/*}}}*/
} Functional;

/*}}}*/

/* Globals */ /*{{{*/
boolean VERBOSE;
boolean MOREVERBOSE;
boolean QUIET;
/* Whether there is need to do any Lua calls during 
 * simulation. */
boolean LUA_CALLS;
/*}}}*/

#endif
