/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_core.c
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

/* Includes */ /*{{{*/
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include "grapham_math.h"
#include "grapham_rand.h"
#include "lua_tools.h"
#include "lua_math.h"
#include "grapham_types.h"
#include "grapham_io.h"

#ifdef _HAVE_NUMLUA
#include <luamatrix.h>
#include <luacomplex.h>
#include <luaspfun.h>
#include <luarng.h>
#include "matrix_type.h"
#endif

#include "custom_type.h"
#include "vector_type.h"
#include "number_type.h"

/*}}}*/

/* Call a density function and return the value */
double call_density(lua_State *L, Node* node) { /*{{{*/
/*
 * In:
 *   L    -- Lua state.
 *   node -- The node whose density is evaluated. Fields used:
 *     value         -- The value of the node.
 *     Nparents      -- Number of parent variables.
 *     is_cdensity   -- Boolean, whether there is a C-density
 *                      to call.
 *     cdensity_args -- Vector of pointers to each argument.
 *     cdensity_lens -- Vector of lengths of each argument.
 *     lua_value     -- The index of the Lua value.
 *     parents       -- The parent nodes. Field used:
 *     censored      -- Whether there is additional censoring.
 *     limits        -- The censoring limits.
 * 
 * Out:
 *   The current value of the log-density of the node.
 */
  int j, Na = (node->Nparents+1);
  double d_value;
  if (node->censored) {
    /* Check first if falls outside censoring values */
    if (node->value[0]<node->limits[0] 
        || node->value[0]>node->limits[1]) 
      return(GRAPHAM_NINF);
  }
  if (node->is_cdensity) {
    /* Call the C-density, if there is one. */
    d_value = (*node->cdensity)(node->cdensity_args, node->cdensity_lens, Na);
  } else {
    /* The function */
    lua_rawgeti(L, LUA_REGISTRYINDEX, node->density);
    /* The arguments */
    if (node->type == NUMBER_TYPE) {
      lua_pushnumber(L, node->value[0]);
    } else {
      lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
    }
    for (j=0; j<node->Nparents; j++) {
      if (node->parents[j]->type == NUMBER_TYPE) {
        lua_pushnumber(L, node->parents[j]->value[0]);
      } else {
        lua_rawgeti(L, LUA_REGISTRYINDEX, node->parents[j]->lua_value);
      }
    }
    lua_call(L, Na, 1);
    /* The return value */
    d_value = lua_tonumber(L,-1);
    lua_pop(L,1); 
  }
  return(d_value);
}

/*}}}*/

/* Compute the initial density values, and check for validity (!=GRAPHAM_NINF) */
void set_initial_likelihoods(lua_State *L, Model* model) { /*{{{*/
/*
 * In:
 *   L     -- Lua state.
 *   model -- The model. Fields used:
 *     N     -- Number of nodes.
 *     nodes -- The nodes. Fields used:
 *       Nparents    -- Number of parent variables.
 *       is_cdensity -- Boolean, whether there is a C-density
 *                      to call.
 *       value       -- The value of the node.
 *       dim         -- Dimension.
 *       parents     -- The parent nodes.
 * 
 * Out:
 *   model -- The updated model. Field updated:
 *     nodes -- The updated nodes. Fields updated:
 *       cdensity_args -- Vector of pointers to each argument.
 *       cdensity_lens -- Vector of lengths of each argument.
 *       density_value -- The initial value of the log-density.
 */
  int j, k, Na;
  double p;
  Node* node;
  for (k=0; k<model->N; k++) {
    node = &(model->nodes[k]);
    if (MOREVERBOSE) fprintf(stderr, " * Likelihood of '%s'\n",node->vname);
    /* Do not evalute, if not necessary: */
    if (node->instantiated && node->Nparents == 0) {
      p = 0;
    } else {
      if (node->is_cdensity) {
        /* Initialise the arguments for the C density function */
        Na = node->Nparents+1;
        node->cdensity_args = (double **)malloc(sizeof(double *)*Na);
        node->cdensity_lens = (int *)malloc(sizeof(int)*Na);
        node->cdensity_args[0] = node->value; 
        node->cdensity_lens[0] = node->dim;
        for (j=0; j<node->Nparents; j++) {
          node->cdensity_args[j+1] = node->parents[j]->value;
          node->cdensity_lens[j+1] = node->parents[j]->dim;
        }
      }
      p = call_density(L, node);
    }
    if (!isfinite(p)) {
      error(L,"The initial log-density of '%s' evaluated -infinity.\n"
            "Please specify valid initial values.\n",node->vname);
    } else {
      node->density_value = p;
    }
  }
}

/*}}}*/

/* Initialise the Lua variables. */
int init_variables(lua_State *L, Model* model) { /*{{{*/
/*
 * In:
 *   L     -- Lua state.
 *   model -- The model. Fields used:
 *     N     -- Number of nodes.
 *     nodes -- The nodes. Fields used:
 *       vname -- Variable name.
 *       type  -- The node type.
 * 
 * Out:
 *   model -- The updated model. Field updated:
 *     nodes -- The updated nodes. Fields updated:
 *       lua_value -- The index to the Lua value.
 * 
 * The function creates a global Lua variable having the same name
 * as in the model, and stores the index to it.
 */
  int k;
  Node* node;
  for (k=0; k<model->N; k++) {
    node = &(model->nodes[k]);
    lua_pushstring(L, node->vname);
    lua_rawget(L, LUA_GLOBALSINDEX);
    if (!lua_isnil(L,-1)) 
      error(L,"Error: '%s': global variable '%s' already defined",
            MODEL_VNAME, node->vname);
    lua_pop(L,1);
    
    switch(node->type) {
    case NUMBER_TYPE:
      number_create(L, node);
      break;
    case VECTOR_TYPE:
      vector_create(L, node);
      break;
#ifdef _HAVE_NUMLUA
    case MATRIX_TYPE:
      matrix_create(L, node);
      break;
#endif
    case CUSTOM_TYPE:
      custom_create(L, node);
      break;
    case UNKNOWN_TYPE:
      ;
    }
    /* Get the created variable */
    lua_getglobal(L, node->vname);
    /* Make a reference to the variable (and pop the value off the stack) */
    node->lua_value = luaL_ref(L, LUA_REGISTRYINDEX);
  }
  return EXIT_SUCCESS;
}

/*}}}*/

/* Call the functional, and increment the mean */
void call_functional(lua_State* L, Functional* functional) { /*{{{*/
/*
 * In:
 *   L          -- Lua state.
 *   functional -- The functional. Fields used:
 *     called     -- The number of times the function is called.
 *     mean       -- The cumulated mean of the functional.
 *     mean_N     -- The length of the mean vector.
 *     is_cfun    -- Whether the functional is in C.
 *     cfun       -- The C function.
 *     cfun_tmp   -- Temporary storage of length mean_N.
 *     cfun_Nargs -- Number of arguments to the function.
 *     cfun_args  -- The list of pointers to the arguments.
 *     cfun_lengths -- The list of lengths of the arguments.
 *     fun        -- The Lua function.
 *     type       -- The type the Lua function returns.
 *
 * Out:
 *   functional -- The functional with updated fields:
 *     mean       -- The mean cumulated by current value.
 *     called     -- Incremented the number functional was called.
 */
  double alpha;
  int i, N;
  long c;
#ifdef _HAVE_NUMLUA
  lua_Matrix* M;
#endif
  
  if (functional->called==0) { /*{{{*/
    if (functional->is_cfun) {
      functional->mean = 
      (double *)malloc(sizeof(double)*functional->mean_N);
      functional->cfun_tmp = 
      (double *)malloc(sizeof(double)*functional->mean_N);
      (functional->cfun)(functional->mean, functional->mean_N,
                          functional->cfun_args, 
                          functional->cfun_lengths,
                          functional->cfun_Nargs);
    } else {
      /* Initialise: call the function */
      lua_rawgeti(L, LUA_REGISTRYINDEX, functional->fun);
      lua_call(L, 0, 1);
      /* Check the output & allocate mem. */
      if (lua_isnil(L, -1)) {
        N = 0;
        functional->type = UNKNOWN_TYPE;
      } else if (lua_isnumber(L, -1)) {
        N = 1;
        functional->type = NUMBER_TYPE;
      } else if (lua_istable(L, -1)) {
        N = lua_objlen(L,-1);
        functional->type = VECTOR_TYPE;
#ifdef _HAVE_NUMLUA
      } else if (lua_isuserdata(L, -1)) {
        functional->type = MATRIX_TYPE;
        M = (lua_Matrix *)lua_touserdata(L, -1);
        if (M->dim == 1) {
          N = M->size;
        } else {
          N = M->size*M->level[0]->size;
        }
#endif
      } else {
        functional->type = CUSTOM_TYPE;
        for (N=1; ; N++) {
          lua_pushnumber(L, N);
          lua_gettable(L, -2);
          if (lua_isnoneornil(L, -1)) {
            lua_pop(L, 1); N--; break;
          }
          lua_pop(L, 1);
        }
      }
      
      functional->mean_N = N;
      functional->mean = (double *)malloc(sizeof(double)*N);
      /* Get the value */
      switch(functional->type) {
      case NUMBER_TYPE:
        functional->mean[0] = luaL_checknumber(L, -1);
        break;
      case VECTOR_TYPE:
        for (i=0; i<N; i++) {
          lua_rawgeti(L, -1, i+1);
          functional->mean[i] = luaL_checknumber(L, -1);
          lua_pop(L,1);
        }
        break;
#ifdef _HAVE_NUMLUA
      case MATRIX_TYPE:
        M = (lua_Matrix*)lua_touserdata(L, -1);
        memcpy(functional->mean, M->data, sizeof(double)*N);
        break;
#endif
      case CUSTOM_TYPE:
        for (i=0; i<N; i++) {
          lua_pushnumber(L, i+1);
          lua_gettable(L, -2);
          functional->mean[i] = luaL_checknumber(L, -1);
          lua_pop(L,1);
        }
        break;
      case UNKNOWN_TYPE:
        break;
      }
      lua_pop(L,1);
    }
    functional->called++;
    if (MOREVERBOSE) 
    fprintf(stderr, " * Called functional the first time.\n");
    /*}}}*/
  } else { /*{{{*/
    c = functional->called;
    N = functional->mean_N;
    alpha =  1.0/((double)c + 1.0);
    functional->called = c+1;
    if (functional->is_cfun) {
      (functional->cfun)(functional->cfun_tmp, N,
                         functional->cfun_args, 
                         functional->cfun_lengths,
                         functional->cfun_Nargs);
      for (i=0; i<N; i++) {
        functional->mean[i] = functional->mean[i]*(1-alpha) 
        + alpha*functional->cfun_tmp[i];
      }
      
    } else {
      lua_rawgeti(L, LUA_REGISTRYINDEX, functional->fun);
      lua_call(L, 0, 1);
      switch(functional->type) {
#ifdef _HAVE_NUMLUA
      case MATRIX_TYPE:
        M = (lua_Matrix *)lua_touserdata(L, -1);
        for (i=0; i<N; i++) {
          functional->mean[i] = functional->mean[i]*(1-alpha) 
          + alpha*M->data[i];
        }
        break;
#endif
      case NUMBER_TYPE:
        functional->mean[0] = functional->mean[0]*(1-alpha) 
        + alpha*lua_tonumber(L, -1);
        break;
      case VECTOR_TYPE:
        for (i=0; i<N; i++) {
          lua_rawgeti(L, -1, i+1);
          functional->mean[i] = functional->mean[i]*(1-alpha) 
          + alpha*lua_tonumber(L, -1);
          lua_pop(L,1);
        }
        break;
      case CUSTOM_TYPE:
        for (i=0; i<N; i++) {
          lua_pushnumber(L, i+1);
          lua_gettable(L, -2);
          functional->mean[i] = functional->mean[i]*(1-alpha) 
          + alpha*lua_tonumber(L, -1);
          lua_pop(L,1);
        }
        break;
      case UNKNOWN_TYPE:
        break;
      }
      lua_pop(L,1);
    }
  } /*}}}*/

}

/*}}}*/

/* Set one Lua global value. */
void set_lua_value(lua_State *L, double value, 
                   Node* node, const int ind) { /*{{{*/
  switch(node->type) {
  case NUMBER_TYPE:
    number_set_single_value(L, value, node);
    break;
  case VECTOR_TYPE:
    vector_set_single_value(L, value, node, ind);
    break;
#ifdef _HAVE_NUMLUA
  case MATRIX_TYPE:
    matrix_set_single_value(L, value, node, ind);
    break;
#endif
  case CUSTOM_TYPE:
    custom_set_single_value(L, value, node, ind);
    break;
  case UNKNOWN_TYPE:
    ;
  }
}

/*}}}*/

/* Update the Lua values of a block. */
void set_node_values(lua_State* L, Block* block) { /*{{{*/
  int k, ind;
  Node* node;
  if (block->Nvars == 1) { /*{{{*/
    /* Just a speedup in this case... */
    node = block->components[0].node;
    for (k=0; k<block->dim; k++) {
      ind = block->components[k].index;
      node->value[ind] = block->value[k];
      /*if (node->kind == SYMMETRIC_KIND) {
        sz = node->size[0];
        i = (ind%sz); j = (ind/sz);
        node->value[i*sz+j] = block->value[k];
      }*/
    }
    if (LUA_CALLS) {
      switch (node->type) {
      case NUMBER_TYPE:
        number_set_single_value(L, block->value[0], node);
        break;
#ifdef _HAVE_NUMLUA
      case MATRIX_TYPE:
        matrix_set_node_block(L, node, block);
        break;
#endif
      case CUSTOM_TYPE:
        custom_set_node_block(L, node, block);
        break;
      case VECTOR_TYPE:
        vector_set_node_block(L, node, block);
        break;
      case UNKNOWN_TYPE:
        break;
      } 
    } /*}}}*/
  } else {
    for (k=0; k<block->dim; k++) {
      node = block->components[k].node;
      ind = block->components[k].index;
      node->value[ind] = block->value[k];
      /*if (node->kind == SYMMETRIC_KIND) {
        sz = node->size[0];
        i = ind%sz; j = ind/sz;
        node->value[i*sz + j] = block->value[k];
      }*/
      if (LUA_CALLS) {
        set_lua_value(L, block->value[k], node, ind);
      }
    }
  }
}

/*}}}*/

/* Draw a proposal and compute the log-density ratio. */
double block_propose(lua_State* L, Block* block, double* chol, double scaling,
                     void (*rand)(double*, const int)) { /*{{{*/
/* In:
 *   L      -- Lua state (passed to set_node_value and call_density).
 *   block  -- The sampling block. Fields used:
 *     dim        -- The dimension.
 *     components -- The components in the block.
 *       index      -- The index of the component.
 *       node       -- The node of the component.
 *         kind          -- The kind of the node.
 *         value         -- The current value of the node.
 *     variables, children -- The variables and the children variables
 *                            of the block.
 *       density_value -- The current value of the log-density.
 *     
 *   chol    -- The upper-triangular scaling matrix.
 *   scaling -- The additional scaling factor.
 *   rand    -- A function generating random variables from the proposal
 *              distribution.
 * 
 * Out:
 *   block -- The updated sampling block. Fields changed:
 *     rand_value     -- The random variable (unscaled) from the proposal.
 *     old_value      -- The old value of the block, before proposing.
 *     proposal_value -- The new proposed value of the block.
 *     variables.new_density_value, 
 *     children.new_density_value   -- The new (proposed) log-density value.
 * 
 *   The return value is the logarithm of the ratio of the densities
 *   of the proposed and the current value, respectively.
 */
  int j, ind; 
  Node* node;
  double old_val, r_val, alpha, p;

  /* Get a random vector from the proposal distribution and scale
   * it with an upper-triangular Matrix and a scalar.*/
  (*rand)(block->rand_value, block->dim);
  mvrand(chol, scaling, block->rand_value, block->proposal_value, block->dim, TRUE);
  
  for (j=0; j<block->dim; j++) {
    node = block->components[j].node;
    ind = block->components[j].index;
    old_val = node->value[ind];
    /* Store the old value, which is reverted if the proposal is not
     * accepted. */
    block->old_value[j] = old_val;
    /* Add the current value of the proposal increment */
    r_val = block->proposal_value[j];
    if (node->kind == INTEGER_KIND) {
      block->proposal_value[j] = old_val + SIGNUM(r_val)*(1.0+floor(fabs(r_val)));
    } else { /* REAL_KIND */
      block->proposal_value[j] = old_val + r_val;
    }
  }
  /* Set the node values of the proposal */
  block->value = block->proposal_value;
  set_node_values(L, block);
  
  alpha = 0.0;
  for (j=0; j<block->Nvars; j++) { 
    node = block->variables[j];
    /* Compute the likelihood difference */
    p = call_density(L,node);
    if (!isfinite(p)) return GRAPHAM_NINF;
    alpha += p - node->density_value;
    node->new_density_value = p;
  }
  for (j=0; j<block->Nchildren; j++) { 
    node = block->children[j];
    /* Compute the likelihood difference */
    p = call_density(L,node);
    if (!isfinite(p)) return GRAPHAM_NINF;
    alpha += p - node->density_value;
    node->new_density_value = p;
  }
  return alpha;
} /*}}}*/

/* Update node information when a proposal is accepted. */
void block_accept(lua_State* L, Block* block) { /*{{{*/
/*
 * In:
 *   L     -- Lua state.
 *   block -- The sampling block. Fields used:
 *     variables -- The variables in the block. Field used:
 *       new_density_value -- The density value with the proposal.
 *     children  -- The children variables of the block. Field used:
 *       new_density_value -- The density value with the proposal.
 * 
 * Out:
 *   block -- The updated sampling block. Fields updated:
 *     accepted  -- The counter of accepted proposals.
 *     variables.density_value,
 *     children.density_value -- The updated density values.
 */
  int j;
  Node* node;
  
  block->accepted++;
  for (j=0; j<block->Nvars; j++) { 
    node = block->variables[j];
    node->density_value = node->new_density_value;
  }
  for (j=0; j<block->Nchildren; j++) { 
    node = block->children[j];
    node->density_value = node->new_density_value;
  }
}
/*}}}*/

/* Update node information when a proposal is rejected. */
void block_reject(lua_State* L, Block* block) { /*{{{*/
/*
 * In:
 *   L     -- Lua state.
 *   block -- The sampling block. Field used:
 *     old_value -- The old value of the block.
 * 
 * Out:
 *   block -- The updated sampling block. Fields updated:
 *     rejected  -- The counter of rejected proposals.
 *     value     -- The pointer to the current value of the block.
 */
  block->rejected++;
  /* Switch value<->old_value */
  block->value = block->old_value;
  set_node_values(L, block);
} /*}}}*/

/* Compute the adaptation weight step size. */
double adapt_step_size(lua_State* L, int k, int fun, boolean is_fun, 
                       double exponent) { /*{{{*/
/*
 * In:
 *   L    -- Lua state.
 *   k    -- The iteration (non-negative).
 *   is_fun   -- Whether there is a custom adaptation weight function.
 *   fun      -- The custom adaptation weight function.
 *   exponent -- The exponent of the adaptation weight.
 * 
 * Out:
 *   The non-negative adaptation weight.
 */
  double eta_k;
  if (is_fun == TRUE) {
    lua_rawgeti(L, LUA_REGISTRYINDEX, fun);
    lua_pushnumber(L, k);
    lua_call(L, 1, 1);
    eta_k = lua_tonumber(L, -1);
    lua_pop(L, 1);
  } else {
    /* The default adaptation sequence */
    if (exponent == 1.0) {
      eta_k = 1.0/((double)k+2.0);
    } else {
      eta_k = 1.0/pow((double)k + 2.0, exponent);
    }
  }
  return eta_k;
} /*}}}*/

/* Standard accept/reject step for a single block. */
double block_accept_reject(lua_State *L, Block* block, Parameters* para,
                           double* chol, double scaling) { /*{{{*/
  double alpha;
  
  alpha = block_propose(L, block, chol, scaling, para->rand);
  
  if (!isfinite(alpha)) {
    alpha = 0;
    goto reject;
  } else {
    /* alpha = min(1,exp(alpha)) */
    alpha = exp(alpha); alpha = (alpha>1)?1:alpha;
  }
  
  /* Do accept-reject */
  if (randu() <= alpha) {
    block_accept(L, block);
  } else {
    reject:
    block_reject(L, block);
  }
  
  return (alpha);
}

/*}}}*/

/* Two-stage delayed rejection step. */
double delayed_accept_reject(lua_State *L, Block* block, Parameters* para, 
                             double *chol, double scaling) { /*{{{*/
  double dr_scaling = para->dr_scaling;
  double alpha, alpha2, p;
  double *u = block->rand_value,
         *u1 = block->tmp_value,
         *y1 = &(block->tmp_value[block->dim]);

  alpha = block_propose(L, block, chol, scaling, para->rand);
  
  if (!isfinite(alpha)) {
    alpha = 0;
    goto reject;
  } else {
    /* alpha = min(1,exp(alpha)) */
    alpha = exp(alpha); alpha = (alpha>1)?1:alpha;
  }
  
  /* Do accept-reject */
  if (randu() <= alpha) {
    block_accept(L, block);
  } else {
    reject:
    /* Save the first proposed value. */
    memcpy(u1, u, sizeof(double)*block->dim);
    if (para->alg == ALG_ASHM) {
      memcpy(y1, block->proposal_value, sizeof(double)*block->dim);
    }
    /* "Reject," i.e. revert the values. */
    block_reject(L, block);
    /* Propose another move with smaller scaling. */
    alpha2 = block_propose(L, block, chol, dr_scaling*scaling, para->rand);
    
    if (!isfinite(alpha2)) {
      alpha2 = 0;
      goto reject2;
    } else {
      /* alpha = min(1,exp(alpha)) */
      alpha2 = exp(alpha2); 
    }
    
    if (alpha2<=alpha) goto reject2;
    /* p = q(y2, u1)/q(x, u1) */
    p = (*para->dram_proposal_ratio)(u1, u, dr_scaling, block->dim);

    if (randu() <= (alpha2-alpha)/(1-alpha)*p) {
      block_accept(L, block);
    } else {
      reject2:
      block_reject(L, block);
    }
    if (para->alg == ALG_ASHM) {
      memcpy(block->proposal_value, y1, sizeof(double)*block->dim);
      memcpy(u, u1, sizeof(double)*block->dim);
    }
    
  }
  
  /* Return the *first-stage acceptance probability!* */
  return (alpha);
}

/*}}}*/

/* Scaling adaptation with the default algorithm */
double scaling_adapt(lua_State* L, double sc, double eta, double alpha, int dim, int k, 
                     Parameters* para) { /*{{{*/
  double acc_opt;
  acc_opt = (dim == 1) ? para->acc_opt1 : para->acc_opt2;
  
/*  return sc*(1+(alpha/acc_opt - 1.0)*eta);*/
  return sc*exp((alpha-acc_opt)*eta);
/*  return sc*(1+(alpha/acc_opt - 1)/sqrt((double)k+2.0));*/
/*  double mu;
  if (alpha<acc_opt) {
    mu = alpha/acc_opt-1.0;
  } else {
    mu = (alpha-acc_opt)/(1.0-acc_opt);
  }
  return sc*(1+mu/((double)k+2.0));*/
}

/*}}}*/

/* Scaling adaptation with a user-defined Lua function */
double scaling_adapt_custom(lua_State* L, double sc, double eta, double alpha, int dim, 
                            int k, Parameters* para) { /*{{{*/
  double sc_;
  lua_rawgeti(L, LUA_REGISTRYINDEX, para->custom_scaling_fun);
  lua_pushnumber(L, sc);
  lua_pushnumber(L, alpha);
  lua_pushnumber(L, dim);
  lua_pushnumber(L, k);
  lua_call(L, 4, 1);
  sc_ = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return sc_;
}

/*}}}*/

/* Update the mean and Cholesky factor of a sampling block */
void block_mean_chol_update(Block* block, Parameters* para, double eta,
                            double alpha, double sc) { /*{{{*/
/* In:
 *   block -- The sampling block to be updated. Fields used:
 *     adapt_chol     -- The Cholesky factor to be updated.
 *     adapt_mean     -- The mean to be updated.
 *     value          -- The new value.
 *     dim            -- The dimension.
 *   para  -- Parameters struct. (not used).
 *   eta   -- Adaptation weight.
 *   alpha -- Observed acceptance probability (not used).
 *   sc    -- The additional scaling factor used when proposing (not used).
 * 
 * Out:
 *   block -- The updated sampling block. Fields changed:
 *     adapt_chol     -- The updated Cholesky factor.
 *     adapt_mean     -- The updated mean.
 *     value          -- Temporary data.
 *     tmp_value      -- Temporary data.
 */
  
  /* Accumulate mean & cholesky factor: */
  sub_vector(block->value, block->value, block->adapt_mean, block->dim);
  /* Update mean */
  mac_vector(block->adapt_mean, block->adapt_mean, eta, block->value, block->dim);
  
  /* Prescale Cholesky factor L down to U so that 
   * U*U' = (1-eta)*L*L' */
  scale_triu(block->adapt_chol, sqrt(1-eta), block->dim);
  /* Update Cholesky factor so that L*L' = U*U' + eta*(x-m)(x-m)' */
  chol_update_linpack(block->adapt_chol, eta, block->value, 
                      block->tmp_value, block->dim);
} /*}}}*/

/* "Rao-Blackwellised" update of the mean and Cholesky factor of a sampling block */
void block_rb_mean_chol_update(Block* block, Parameters* para, 
                               double eta, double alpha, double sc) { /*{{{*/
/* In:
 *   block -- The sampling block to be updated. Fields used:
 *     adapt_chol     -- The Cholesky factor to be updated.
 *     adapt_mean     -- The mean to be updated.
 *     proposal_value -- The proposed value.
 *     old_value      -- The old value.
 *     dim            -- The dimension.
 *   para  -- Parameters struct. (not used).
 *   eta   -- Adaptation weight.
 *   alpha -- Observed acceptance probability (first-stage, if delayed
 *            rejection is used.)
 *   sc    -- The additional scaling factor used when proposing.
 * 
 * Out:
 *   block -- The updated sampling block. Fields changed:
 *     adapt_chol     -- The updated Cholesky factor.
 *     adapt_mean     -- The updated mean.
 *     proposal_value -- Temporary data.
 *     old_value      -- Temporary data.
 *     tmp_value      -- Temporary data.
 */
  double alpha1;
  
/*  eta = adapt_step_size(L, k, para);*/
  alpha1 = 1-alpha;
  
  /* Accumulate mean & cholesky factor: */
  sub_vector(block->proposal_value, 
             block->proposal_value, block->adapt_mean, block->dim);
  sub_vector(block->old_value,
             block->old_value,      block->adapt_mean, block->dim);
  /* Update mean */
  mac_vector(block->adapt_mean, 
             block->adapt_mean, alpha*eta,  block->proposal_value, block->dim);
  mac_vector(block->adapt_mean,
             block->adapt_mean, alpha1*eta, block->old_value, block->dim);
  
  /* Prescale Cholesky factor L down to U so that 
   * U*U' = (1-eta)*L*L' */
  scale_triu(block->adapt_chol, sqrt(1-eta), block->dim);
  /* Update Cholesky factor so that L*L' = U*U' + eta*(x-m)(x-m)' */
  chol_update_linpack(block->adapt_chol, alpha*eta, block->proposal_value, 
                      block->tmp_value, block->dim);
  chol_update_linpack(block->adapt_chol, alpha1*eta, block->old_value, 
                      block->tmp_value, block->dim);
}
/*}}}*/

/* Adaptive direction scaling update of the Cholesky factor of a sampling block. */
void block_ashm_chol_update(Block* block, Parameters* para, 
                            double eta, double alpha, double sc) { /*{{{*/
/* In:
 *   block -- The sampling block to be updated. Fields used:
 *     adapt_chol     -- The Cholesky factor to be updated.
 *     proposal_value -- The proposed value.
 *     old_value      -- The old value.
 *     rand_value     -- The unscaled proposal increment.
 *     dim            -- The dimension.
 *   para  -- Parameters struct. Fields used:
 *     acc_opt1, acc_opt2 -- Desired acceptance probabilities.
 *   eta   -- The adaptation weight.
 *   alpha -- Observed acceptance probability (first-stage, if delayed
 *            rejection is used.)
 *   sc    -- The additional scaling factor used when proposing.
 * 
 * Out:
 *   block -- The updated sampling block. Fields changed:
 *     adapt_chol     -- The updated Cholesky factor.
 *     proposal_value -- Temporary data.
 *     tmp_value      -- Temporary data.
 */

  int i, info;
  double s;
  double *u = block->rand_value;
  double acc_opt = (block->dim == 1) ? para->acc_opt1 : para->acc_opt2;
  double a = alpha-acc_opt;
  double eta_d = MIN(eta * block->dim, 1);

  /* Compute the unscaled proposal vector length */
  s = 0.0;
  for (i=0; i<block->dim; i++) s += u[i]*u[i];
  s = (s > 0) ? s : 1.0;
  s *= sc*sc;
  sub_vector(block->proposal_value, block->proposal_value, block->old_value, block->dim);
  
  /* Update Cholesky factor so that L*L' = U*U' + eta*a*u*u' */
  if (a>=0) {
    chol_update_linpack(block->adapt_chol, a*eta_d/s, block->proposal_value, 
                        block->tmp_value, block->dim);
  } else { /* a<0 */
    info = chol_downdate_linpack(block->adapt_chol, -a*eta_d/s, 
                                 block->proposal_value, 
                                 block->tmp_value, block->dim);
    if (VERBOSE && info != 0) 
      fprintf(stderr, "Warning: Cholesky factor downdate failed!\n");
  }
}
/*}}}*/

/* One full MCMC step; update all blocks sequentially, and perform
 * the adaptations. */
void mcmc_step(lua_State *L, Model* model, Parameters* para,
                     const int k) { /*{{{*/
  int n, j;
  Block *block, *block_;
  double eta = 0.5, eta_sc = 0.5, alpha, scaling, *chol, p_mix = 0.0;
  
  /* Compute step sizes */
  if (para->mean_chol_update != NULL) {
    eta = adapt_step_size(L, k, para->adapt_weight_fun, 
                          para->is_adapt_weight_fun, 
                          para->adapt_weight_exp);
  }
  if (para->scaling_adapt != NULL) {
    eta_sc = adapt_step_size(L, k, para->adapt_weight_sc_fun, 
                             para->is_adapt_weight_sc_fun, 
                             para->adapt_weight_sc_exp);
  }
  
  n = model->Nblocks;
  while (n >= 1) {
    if (para->random_scan) {
      /* Fisher-Yates shuffle */
      j = randint(n-1);
      n--;
      block_ = model->blocks[n];
      block = model->blocks[j];
      model->blocks[n] = block;
      model->blocks[j] = block_;
    } else {
      block = model->blocks[--n];
    }
    
    /* If the mixed proposal is used */
    if (para->is_mix_weight_fun) {
      lua_rawgeti(L, LUA_REGISTRYINDEX, para->mix_weight_fun);
      lua_pushnumber(L, k);
      lua_call(L, 1, 1);
      p_mix = lua_tonumber(L, -1);
      lua_pop(L, 1);
    } 
    
    /* Do the accept-reject step */
    if ( (para->init == INIT_TRAD && k<para->nburn) 
         || (para->is_mix_weight_fun && randu() < p_mix) ) {
      chol = block->init_chol;
      scaling = block->init_scaling;
    } else {
      chol = block->adapt_chol;
      scaling = block->scaling;
    }
    alpha = (*para->accept_reject)(L, block, para, chol, scaling);

    /* Do adaptations */
    if (!para->init == INIT_FREEZE || k<para->nburn) {
      if (para->scaling_adapt != NULL) {
        block->scaling = 
        (*para->scaling_adapt)(L, block->scaling, eta_sc, alpha, block->dim, k, para);
      }
      if (para->mean_chol_update != NULL) {
        (*para->mean_chol_update)(block, para, eta, alpha, scaling);
      }
    }
  }
  
} /*}}}*/

/* Initial setup for different adaptation & sampling strategies. */
void activate_scaling(Parameters* para) { /*{{{*/
    if (para->is_custom_scaling_fun) {
      para->scaling_adapt = *scaling_adapt_custom;
    } else {
      para->scaling_adapt = *scaling_adapt;
    }  
}
/*}}}*/

/* Set the adaptation functions corresponding to a selected algorithm,
 * and do the blocking of variables (that are not beforehand assigned
 * to a block).*/
void algorithm_setup(lua_State *L, Model* model, Parameters* para) { /*{{{*/
  int k;
  
  /* Setup of different algorithms */
  switch(para->alg) {
  case ALG_ASCM:
    para->mean_chol_update = NULL;
    activate_scaling(para);
    break;
  case ALG_AM:
    para->mean_chol_update = *block_mean_chol_update;
    para->scaling_adapt = NULL;
    break;
  case ALG_AMS:
    para->mean_chol_update = *block_mean_chol_update;
    activate_scaling(para);
    break;
  case ALG_RBAM:
    para->mean_chol_update = *block_rb_mean_chol_update;
    para->scaling_adapt = NULL;
    break;
  case ALG_RBAMS:
    para->mean_chol_update = *block_rb_mean_chol_update;
    activate_scaling(para);
    break;
  case ALG_ASHM:
    para->mean_chol_update = *block_ashm_chol_update;
    para->scaling_adapt = NULL;
    break;
  case ALG_METROPOLIS:
    para->mean_chol_update = NULL;
    para->scaling_adapt = NULL;
    break;
  }
  
  /* Whether to activate delayed rejection */
  if (para->dr_scaling>0.0) {
    for (k=0; k<model->N; k++) 
      if (model->nodes[k].kind == INTEGER_KIND)
        error(L, "Error: Delayed Rejection not applicable with integer variables.");
    para->accept_reject = *delayed_accept_reject;
  } else {
    para->accept_reject = *block_accept_reject;
  }
  
  /* Fill the remaining blocks as requested */
  switch(para->blocking) {
  case BLOCKING_SC:
    fill_single_components(L, model);
    break;
  case BLOCKING_NODE:
    fill_blocks(L, model);
    break;
  case BLOCKING_FULL:
    fill_single_block(L, model);
    break;
  }
  
  /* Read and prepare blocks */
  read_blocks(L, model);

}

/*}}}*/

void grapham_preinit(lua_State *L) { /*{{{*/
  FILE* tmpf;
  
  /* Open Lua and load standard libraries */
  if (VERBOSE) fprintf(stderr, "Opening Lua standard libraries...\n");
  luaL_openlibs(L);

  /* Load init script */
  if (VERBOSE) fprintf(stderr, "Loading the Lua init script...\n");
#ifdef _GRAPHAM_ROOT
  tmpf = fopen(QUOTE(_GRAPHAM_ROOT) "/" LUA_INIT_SCRIPT, "r");
  if (tmpf != NULL) {
    fclose(tmpf);
    if (luaL_loadfile(L, QUOTE(_GRAPHAM_ROOT) "/" LUA_INIT_SCRIPT) 
        || lua_pcall(L, 0, 0, 0))
      error(L, "Error: %s", lua_tostring(L, -1));
  } else 
#endif
  {
    if (luaL_loadfile(L, LUA_INIT_SCRIPT) || lua_pcall(L, 0, 0, 0))
    error(L, "Error: %s", lua_tostring(L, -1));
  }

#ifdef _HAVE_NUMLUA
  if (VERBOSE) fprintf(stderr, "Opening Numlua libraries...\n");
  lua_pushcfunction(L, luaopen_luaspfun);
  lua_pushstring(L, "spfun");
  lua_call(L, 1, 0);
  lua_pushcfunction(L, luaopen_luacomplex);
  lua_pushstring(L, "complex");
  lua_call(L, 1, 0);
  lua_pushcfunction(L, luaopen_luamatrix);
  lua_pushstring(L, "matrix");
  lua_call(L, 1, 0);
  lua_pushcfunction(L, luaopen_luarng);
  lua_pushstring(L, "rng");
  lua_call(L, 1, 0);
#endif 

  if (VERBOSE) fprintf(stderr, "Exporting math functions...\n");
  if (export_lua_math(L) < 0) error(L, "Error: %s", lua_tostring(L, -1));
} /*}}}*/

void grapham_init(lua_State* L, Model *model, Functional* functional, 
                  Parameters* para) { /*{{{*/
  FILE* tmpf;
  
  /* Load the preprocessing script */
  if (VERBOSE) fprintf(stderr, "Loading the Lua preprocessing script...\n");
#ifdef _GRAPHAM_ROOT
  tmpf = fopen(QUOTE(_GRAPHAM_ROOT) "/" LUA_PREPROCESS_SCRIPT, "r");
  if (tmpf != NULL) {
    fclose(tmpf);
    if (luaL_loadfile(L, QUOTE(_GRAPHAM_ROOT) "/" LUA_PREPROCESS_SCRIPT) 
        || lua_pcall(L, 0, 0, 0))
      error(L, "Error: %s", lua_tostring(L, -1));
  } else 
#endif
  {
    if (luaL_loadfile(L, LUA_PREPROCESS_SCRIPT) || lua_pcall(L, 0, 0, 0))
    error(L, "Error: %s", lua_tostring(L, -1));
  }
  
  /* Read the parameters struct */
  if (VERBOSE) fprintf(stderr, "Reading '%s'...\n", PARA_VNAME);
  read_parameters(L, model, para);
  seed_random_generator(para->seed);  

  /* Read the model structure */
  if (VERBOSE) fprintf(stderr, "Reading '%s'...\n", MODEL_VNAME);
  read_model(L, model, para);  

  /* Initialise the variables in Lua */
  if (VERBOSE) fprintf(stderr, "Initialising variables in Lua...\n");
  init_variables(L, model);
  
  if (VERBOSE) fprintf(stderr, "Reading the functional...\n");
  read_functional(L, model, para, functional);
  
  /* Do adaptive algorithm setup (if any) */
  if (VERBOSE) fprintf(stderr, "Setting up the adaptation algorithm...\n");
  algorithm_setup(L, model, para);
  
  /* Compute initial likelihoods */
  if (VERBOSE) fprintf(stderr, "Computing initial likelihoods...\n");
  set_initial_likelihoods(L, model);
  
  /* Check connectivity  acyclicity */
  if (VERBOSE) fprintf(stderr, "Checking model acyclicity  connectivity...\n");
  if (check_dag(model) != EXIT_SUCCESS)
    error(L, "Error: there is a cycle in the model.");
  check_connectivity(model);
  
  if (para->output.write_file) {
    if (VERBOSE) fprintf(stderr, 
                         "Writing output file header...\n");
    prepare_outfile(L, model, &(para->output));
  }

  if (para->adapt_output.write_file) {
    if (VERBOSE) fprintf(stderr, 
                         "Writing adaptation parameter output file header...\n");
    prepare_outadapt(L, model, para, &(para->adapt_output));
  }

} /*}}}*/

void grapham_run(lua_State *L, Model *model, Functional* functional,
                         Parameters* para) { /*{{{*/
  unsigned long k, 
  niter = para->niter,
  nburn = para->nburn, 
  nthin = para->nthin;

  if (VERBOSE) fprintf(stderr, "Starting sampling...\n");
  for (k=0; k<niter+nburn; k++) {
    
    /* Do propose-accept-reject step */
    mcmc_step(L, model, para, k);

    if ((k%nthin) == 0) {
      if (para->adapt_output.write_file) 
        (*para->adapt_output.add_record_outfile)(&(para->adapt_output));
      if (k>=nburn) {
        if (functional->evaluate) call_functional(L, functional);
        /* Write output, if requested */
        if (para->output.write_file) 
          (*para->output.add_record_outfile)(&(para->output));
      }
    }
  }
} /*}}}*/

void grapham_close(lua_State* L, Model *model, Parameters* para) { /*{{{*/
  
  if (para->output.write_file) {
    if (VERBOSE) fprintf(stderr, "Flushing output file...\n");
    fclose(para->output.file);
  }

  if (para->adapt_output.write_file) {
    if (VERBOSE) fprintf(stderr, "Flushing adaptation output file...\n");
    fclose(para->adapt_output.file);
  }

  if (para->outcfg != NULL) {
    if (VERBOSE) fprintf(stderr, "Writing config...\n");
    write_outcfg(para->outcfg, model);
  }
  
  if (para->is_close_fun) {
    if (VERBOSE) fprintf(stderr, "Calling '%s.%s'...\n",
                         PARA_VNAME, PARA_CLOSE_VNAME);
    lua_rawgeti(L, LUA_REGISTRYINDEX, para->close_fun);
    lua_call(L, 0, 0);
  }
} /*}}}*/

