/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_io.c
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
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <dlfcn.h>  /* dlopen, dlsym */
#include <unistd.h> /* getcwd */

#include "grapham_io.h"
#include "lua_tools.h"
#include "grapham_math.h"
#include "grapham_rand.h"

#ifdef _HAVE_NUMLUA
#include "matrix_type.h"
#endif

#include "custom_type.h"
#include "vector_type.h"
#include "number_type.h"

#include "grapham_lib.h"
/*}}}*/

/* Find a given variable name within the model. */
int find_variable(Model* model, const char *str) { /*{{{*/
  int i;
  for (i=0; i<model->N; i++) {
    if (!strcmp(str,model->nodes[i].vname)) {
      return(i);
    }
  }
  return(-1);
}

/*}}}*/

/* Given a string, finds a variable whose name matches,
 * and an index of the vector, if given. (like 'x[1]')
 * On error, the return value is negative. */
int find_variable_index(Model* model, const char *str, 
                  int *i) { /*{{{*/
  char *obrk, *cbrk = NULL, *obrk2 = NULL, *cbrk2 = NULL, *endptr;
  Node* node;
  int k, j, len, ind, str_len = strlen(str);
  
  /* Find the possible first index */
  obrk = strchr(str, '[');
  if (obrk != NULL) {
    len = obrk-str;
    cbrk = strchr(obrk, ']');
    /* Check that closing bracket is found */
    if (cbrk == NULL) return -1;
  } else {
    len = strlen(str);
  }
  
  /* Find the matching variable name */
  ind = -1;
  for (k=0; k<model->N; k++) {
    node = &(model->nodes[k]);
    if (strlen(node->vname) == len 
        && strncmp(str, node->vname, len) == 0) {
      ind = k;
      break;
    }
  }
  
  /* Parse the index (or indices) */
  if (ind>=0) {
    if (obrk != NULL) {
      endptr = cbrk-1;
      *i = strtol(obrk+1, &endptr, 0)-1;
      if (endptr == obrk+1 || *i < 0 || *i > node->dim-1) return -1;
      /* Check if the closing bracket was not the last one */
      if (cbrk != str+str_len-1) {
        /* Then, see if the r.v. is a matrix */
        if (node->Ndims < 2) return -1;
        obrk2 = strchr(cbrk, '[');
        if (obrk2 == NULL) return -1;
        cbrk2 = strchr(obrk2, ']');
        if (cbrk2 == NULL || cbrk2 != str+str_len-1) return -1;
        endptr = cbrk-1;
        j = strtol(obrk2+1, &endptr, 0)-1;
        /* Check the sanity of indices */
        if (endptr == obrk2+1 || j < 0 || j > node->size[1]-1 
            || *i > node->size[0]-1) return -1;
        if (node->kind == SYMMETRIC_KIND && j>*i) {
          *i = node->size[0]*(*i) + j;
        } else {
          *i = node->size[0]*j + (*i);
        }
      }
    } else {
      *i = -1;
    }
  } 
  return ind;
}
/*}}}*/

/* Open output file for writing. */
FILE* open_outfile(lua_State *L, char const *fname, 
                   char const *mode) { /*{{{*/
  FILE* outFile = fopen(fname, mode);
  if (outFile == NULL) {
    error(L, "Error: Cannot open '%s' for writing.", fname);
  } 
  return outFile;
}

/*}}}*/

/* Print the name of one component into output. */
void print_component_name(FILE* ofile, Node* node, int ind, 
                          const char* prefix, const char* postfix) { /*{{{*/
  int sz;
  const char *pref = (prefix == NULL) ? "" : prefix;
  const char *postf = (postfix == NULL) ? "" : postfix;
  if (node->type == NUMBER_TYPE) {
    fprintf(ofile, "%s%s%s", pref, node->vname, postf);
  } else {
    if (node->Ndims == 1) {
      fprintf(ofile, "%s%s[%i]%s", pref, node->vname,  ind+1, postf);
    } else {
      sz = node->size[0];
      fprintf(ofile, "%s%s[%i][%i]%s", pref, node->vname, 
              (ind%sz)+1, (ind/sz)+1, postf);
    }
  }
} /*}}}*/

/* Write the header line with variable names and prepare output. */
void prepare_outfile(lua_State* L, Model* model, 
                     Output* output) { /*{{{*/
  int k,j,i,ind,n;
  const char* str;
  Node* node;
  boolean first=TRUE;
  FILE* ofile = output->file;
  
  lua_getglobal(L, PARA_VNAME);
  lua_getfield(L, -1, PARA_OUTVARS_VNAME);
  if (lua_isnil(L, -1)) {
    /* The default: save all uninstantiated variables */
    lua_pop(L, 1);
    lua_newtable(L);
    n = 0;
    for (k=0;k<model->N; k++) {
      node = &(model->nodes[k]);
      if (node->instantiated) continue;
      lua_pushstring(L, node->vname);
      lua_rawseti(L, -2, ++n);
    }
    lua_setfield(L, -2, PARA_OUTVARS_VNAME);
    lua_getfield(L, -1, PARA_OUTVARS_VNAME);
  } else if (!lua_istable(L, -1)) {
    error(L, "Error: '%s.%s' should be a table.", PARA_VNAME, PARA_OUTVARS_VNAME);
  }

  output->vars = (double **)malloc(sizeof(double *)*model->dim);
  n = 0;
  for (k=0; k<lua_objlen(L, -1); k++) {
    lua_rawgeti(L, -1, k+1);
    if (!lua_isstring(L, -1))
      error(L, "Error: '%s.%s' should be a table of variable names.", 
            PARA_VNAME, PARA_OUTVARS_VNAME);
    str = lua_tostring(L, -1);
    ind = find_variable_index(model, str, &i);
    lua_pop(L, 1); /* str */

    if (ind < 0) 
      error(L, "Error: '%s.%s': Cannot find variable '%s'", 
            PARA_VNAME, PARA_OUTVARS_VNAME, str);
    node = &(model->nodes[ind]);
    if (i < 0) {
      /* No index specified: use all */
      for (j=0; j<node->dim; j++) {
        /* Omit upper-diagonal part for symmetric matrices. */
        if (node->kind == SYMMETRIC_KIND 
            && (j/node->size[0]) > (j%node->size[0])) continue;
        output->vars[n] = &(node->value[j]);
        if (first) {
          first = FALSE;
        } else {
          fprintf(ofile, OUTFILE_FIELD_SEP);
        }
        print_component_name(ofile, node, j, "\"", "\"");
        n++;
      }
    } else {
      if (i >= node->dim) 
        error(L, "Error: '%s.%s': index exceeds variable dimension: '%s'", 
              PARA_VNAME, PARA_OUTVARS_VNAME, str);
      output->vars[n] = &(node->value[i]);
      if (first) {
        first = FALSE;
      } else {
        fprintf(ofile, OUTFILE_FIELD_SEP);
      }
      print_component_name(ofile, node, i, "\"", "\"");
      n++;
    }
  }
  output->Nvars = n;
  output->vars = (double **)realloc(output->vars, sizeof(double *)*n);
  fprintf(ofile,OUTFILE_RECORD_SEP);
}
/*}}}*/

/* Write the header line with adaptation variable names and prepare output. */
void prepare_outadapt(lua_State* L, Model* model, Parameters* para,
                     Output* output) { /*{{{*/
  int k, i, j, N, dim, dim2;
  boolean save_sc = FALSE, save_chol = FALSE, save_mean = FALSE;
  FILE* ofile = output->file;
  Block* block;
  
  switch(para->alg) {
  case ALG_ASCM:
    save_sc = TRUE; 
    break;
  case ALG_AM:
  case ALG_RBAM:
    save_chol = TRUE; save_mean = TRUE;
    break;
  case ALG_AMS:
  case ALG_RBAMS:
    save_sc = TRUE; save_chol = TRUE; save_mean = TRUE;
    break;
  case ALG_ASHM:
    save_chol = TRUE;
    break;
  case ALG_METROPOLIS:
    break;
  }
  dim = 0; dim2 = 0;
  for (k=0; k<model->Nblocks; k++) {
    block = model->blocks[k];
    dim += block->dim;
    /* The elements in the upper-diagonal matrix. */
    dim2 += (block->dim*block->dim + block->dim)/2;
  }
  output->Nvars = 0;
  output->Nvars += save_sc ? model->Nblocks : 0;
  output->Nvars += save_mean ? dim : 0;
  output->Nvars += save_chol ? dim2 : 0;
  
  output->vars = (double **)malloc(sizeof(double *)*output->Nvars);
  N = 0;
  for (k=0; k<model->Nblocks; k++) {
    block = model->blocks[k];
    if (save_sc) {
      if (N>0) fprintf(ofile, ",");
      output->vars[N++] = &(block->scaling);
      fprintf(ofile, "\"SC_");
      for (i=0; i<block->dim; i++) {
        print_component_name(ofile, block->components[i].node,
                             block->components[i].index, (i>0)?"#":NULL, NULL);
      }
      fprintf(ofile, "\"");
    }
    if (save_mean) {
      for (i=0; i<block->dim; i++) {
        if (N>0) fprintf(ofile, ",");
        output->vars[N++] = &(block->adapt_mean[i]);
        print_component_name(ofile, block->components[i].node,
                             block->components[i].index, "\"MEAN_", "\"");
      }
    }
    if (save_chol) {
      for (i=0; i<block->dim; i++) {
        for (j=i; j<block->dim; j++) {
          if (N>0) fprintf(ofile, ",");
          output->vars[N++] = &(block->adapt_chol[IND(i,j,block->dim)]);
          print_component_name(ofile, block->components[i].node,
                               block->components[i].index, "\"CHOL_", NULL);
          print_component_name(ofile, block->components[j].node,
                               block->components[j].index, "#", "\"");
        }
      }
    }
  }
  fprintf(ofile, "\n");
}
/*}}}*/

/* Write one record (line) to the output file. */
void add_ascii_record(Output* output) { /*{{{*/
  int k, ret;
  boolean first=TRUE;
  FILE* ofile = output->file;
  for (k=0;k<output->Nvars; k++) {
    if (first) {
      first = FALSE;
    } else {
      fprintf(ofile, OUTFILE_FIELD_SEP);
    }
    fprintf(ofile, OUTFILE_NUMFORMAT, *(output->vars[k]));
  }
  ret = fprintf(ofile, OUTFILE_RECORD_SEP);
  if (ret < 0) {
      fprintf(stderr, "Error writing output file.\n");
      exit(EXIT_FAILURE);
    }
}

/*}}}*/

/* Write one record (line) to the output file. */
void add_binary_record(Output* output) { /*{{{*/
  int k;
  size_t sz;
  FILE* ofile = output->file;
  for (k=0;k<output->Nvars; k++) {
    sz = fwrite(output->vars[k], sizeof(double), 1, ofile);
    if (sz != 1) {
      fprintf(stderr, "Error writing output file.\n");
      exit(EXIT_FAILURE);
    }
  }
}

/*}}}*/

/* Read a numeric field from a table (if exists) */
double read_numeric_field(lua_State *L, const char* sname, 
                          const char* fname, double default_value) { /*{{{*/
  double value;
  lua_getglobal(L, sname);
  lua_getfield(L, -1, fname);
  if (!lua_isnil(L, -1)) {
    if (!lua_isnumber(L, -1))
      error(L, "Error: '%s.%s' should be a number.",sname,fname);
    value = lua_tonumber(L, -1);
  } else {
    value = default_value;
  }
  lua_pop(L, 2);
  return value;
}

/*}}}*/

/* Check that the node has symmetric value. Return non-zero, if 
 * non-symmetric. In addition, mark lower-diagonal components as
 * if they would be "in block" */
int check_symmetric(Node* node) { /*{{{*/
  int i,j,n,m;
  n = node->size[0]; m = node->size[1];
  if (n != m) return EXIT_FAILURE;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      if (node->value[IND(i,j,n)] != node->value[IND(j,i,n)]) {
        return EXIT_FAILURE;
      }
      node->in_block[IND(i,j,n)] = TRUE;
    }
  }
  return EXIT_SUCCESS;
}
/*}}}*/

/* Try to read value for the node. */
int read_node_value(lua_State *L, Node* node) { /*{{{*/
  if (lua_isnumber(L,-1)) {
    if (node->type != NUMBER_TYPE) return EXIT_FAILURE;
    if (node->kind == INTEGER_KIND) {
      node->value[0] = luaL_checkinteger(L,-1);
    } else {
      node->value[0] = luaL_checknumber(L,-1);
    }
  } else if (node->type == VECTOR_TYPE 
#ifdef _HAVE_NUMLUA
             || node->type == MATRIX_TYPE
#endif
             || node->type == NUMBER_TYPE) {
    if (vector_read_node(L, node) != EXIT_SUCCESS
#ifdef _HAVE_NUMLUA
        && matrix_read_node(L, node) != EXIT_SUCCESS
#endif        
        && custom_read_node(L, node) != EXIT_SUCCESS)
    return EXIT_FAILURE;
  } else {
    if (custom_read_node(L, node) != EXIT_SUCCESS) return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
/*}}}*/

/* Read the model structure */
void read_model(lua_State *L, Model* model, Parameters* para) { /*{{{*/
  int k, i, j, l;
  const char* str;
  size_t str_len;
  Node* node, *parent;

  model->N = 0;
  model->dim = 0;
  model->nodes = NULL;

  k = model->N;

  lua_getglobal(L, MODEL_VNAME);
  if (!lua_istable(L, -1))
    error(L, "Error: '%s' should be a table", MODEL_VNAME);
  
  /* 1) Count the variables in model */
  lua_pushnil(L);
  while (lua_next(L, -2) != 0) {
    model->N++;
    /* removes 'value'; keeps 'key' for next iteration */
    lua_pop(L, 1);
  }
  
  /*printf("%i\n",model -> N);*/
  model->nodes = (Node *)realloc(model->nodes, sizeof(Node)*model->N);
  
  if (VERBOSE) fprintf(stderr, " * Going through %i variables\n",
                       model->N);
  /* 2) Copy the variable names (and check them), and the functions */
  lua_pushnil(L);
  while (lua_next(L, -2) != 0) {
    if (!lua_istable(L, -1) || !lua_isstring(L, -2))
      error(L, "Error: '%s' should consist of (string,table) pairs.",
            MODEL_VNAME);
    /* uses 'key' (at index -2) and 'value' (at index -1) */
    str = lua_tolstring(L, -2, &str_len);
    
    node = &(model->nodes[k]);
    /* Set default values: */
    node->Nchildren = 0;
    /* By default: size[0] == dim */
    node->size = &(node->dim);
    
    /* the variable name */
    node->vname = (char *)malloc(sizeof(char)*str_len+1);
    (void) strcpy(node->vname, str);
  
    /* .instantiated */ 
    if (MOREVERBOSE)
    fprintf(stderr, " --> '%s': %s\n", node->vname, MODEL_INSTANTIATED_VNAME);
    lua_getfield(L, -1, MODEL_INSTANTIATED_VNAME);
    node->instantiated = luaL_checkinteger(L, -1);
    if (node->instantiated != FALSE && node->instantiated != TRUE)
    error(L, "'%s.%s.%s' must be 0 or 1.", MODEL_VNAME, node->vname, 
          MODEL_INSTANTIATED_VNAME);
    lua_pop(L, 1);
    
    /* .density */ /*{{{*/
    if (MOREVERBOSE) 
    fprintf(stderr, " --> '%s': %s\n", node->vname, MODEL_DENSITY_VNAME);
    lua_getfield(L, -1, MODEL_DENSITY_VNAME);
    if (lua_isfunction(L, -1)) {
      node->density = luaL_ref(L, LUA_REGISTRYINDEX);
      node->is_cdensity = FALSE;
      LUA_CALLS = TRUE;
      
    } else if (lua_isstring(L, -1)) {
      str = lua_tostring(L, -1);
      node->cdensity = NULL;
      if (para->clib != NULL) {
        /* Try to find from user-specified library */
        *(void **)(&node->cdensity) = dlsym(para->clib, str);
      } 
      if (node->cdensity == NULL) {
        /* Try to find from grapham_dlib symbols */
        for (j=0; ; j++) {
          if (grapham_dlib[j].name == NULL) {
            break;
          } else {
            if (strcmp(grapham_dlib[j].name, str) == 0) {
              node->cdensity = grapham_dlib[j].c_fun;
            }
          }
        }
      }
      if (node->cdensity == NULL) 
      error(L, "'%s.%s.%s': cannot find symbol '%s'.",
            MODEL_VNAME, node->vname, MODEL_DENSITY_VNAME, str);
      
      node->is_cdensity = TRUE;
      lua_pop(L, 1);
    } else {
      error(L, "Error: '%s.%s.%s' should be a function or a string.", MODEL_VNAME,
            node->vname, MODEL_DENSITY_VNAME);
    } /*}}}*/
    
    /* .type */ /*{{{*/
    if (MOREVERBOSE) fprintf(stderr, " --> '%s': %s\n", node->vname, MODEL_TYPE_VNAME);
    lua_getfield(L, -1, MODEL_TYPE_VNAME);
    if (!lua_isstring(L,-1)) error(L, "'%s.%s' should be a string.", 
                                   MODEL_VNAME, MODEL_TYPE_VNAME);
    str = lua_tostring(L, -1);
    if (!strcmp(str,TYPE_NUMBER_NAME)) node->type = NUMBER_TYPE;
    else if (!strcmp(str,TYPE_VECTOR_NAME)) node->type = VECTOR_TYPE;
    else if (!strcmp(str,TYPE_CUSTOM_NAME)) node->type = CUSTOM_TYPE;
#ifdef _HAVE_NUMLUA
    else if (!strcmp(str,TYPE_MATRIX_NAME)) node->type = MATRIX_TYPE;
#endif
    else error(L, "'%s.%s': Invalid value: '%s'.",
               MODEL_VNAME, MODEL_TYPE_VNAME, str); 
    lua_pop(L, 1); 
    /*}}}*/
    
    /* .dim */ /*{{{*/
    lua_pushstring(L, MODEL_DIM_VNAME);
    lua_rawget(L, -2);
    if (!lua_isnil(L,-1)) {
      if (lua_isnumber(L,-1)) {
        node->Ndims = 1;
        node->dim = (unsigned)luaL_checkinteger(L,-1);
        if (node->type == NUMBER_TYPE && node->dim != 1)
          error(L, "Error: Variable '%s' has type '%s' but '%s'=1.",
                node->vname, TYPE_NUMBER_NAME, MODEL_DIM_VNAME);
      } else if (lua_istable(L,-1)) {
        node->Ndims = lua_objlen(L, -1);
        if (node->Ndims<1 || node->Ndims>2) 
          error(L, "Error: Only variables with 1 or 2 dims are supported.");
        node->size = (int *)malloc(sizeof(int)*node->Ndims);
        node->dim = 1;
        for (j=0; j<node->Ndims; j++) {
          lua_rawgeti(L, -1, j+1);
          node->size[j] = luaL_checkinteger(L, -1);
          node->dim *= node->size[j];
          lua_pop(L, 1);
        }
      } else {
        error(L, "Error: '%s.%s' should be a table or a number.", 
              MODEL_VNAME, MODEL_DIM_VNAME);
      }
    } else {
      if (node->type == CUSTOM_TYPE
#ifdef _HAVE_NUMLUA
          || node->type == MATRIX_TYPE
#endif
          ) 
        error(L, "Error: The type of '%s' requires explicitly "
              "defined dimension.", node->vname);
      node->Ndims = 1;
      node->dim = 1;
    }
    
    lua_pop(L, 1);
    model->dim += node->dim;
    if (MOREVERBOSE) {
      fprintf(stderr, " --> '%s': %s: %i", 
              node->vname, MODEL_DIM_VNAME, node->dim);
      if (node->Ndims>1) {
        fprintf(stderr, " [");
        for (j=0; j<node->Ndims; j++)
          fprintf(stderr, "%i,", node->size[j]);
        fprintf(stderr, "]");
      }
      fprintf(stderr, "\n");
    }
    /*}}}*/

    /* .kind */ /*{{{*/
    if (MOREVERBOSE) fprintf(stderr, " --> '%s': %s\n", 
                             node->vname, MODEL_KIND_VNAME);
    lua_getfield(L, -1, MODEL_KIND_VNAME);
    if (!lua_isstring(L,-1)) 
    error(L, "'%s.%s.%s' should be a string.", 
          MODEL_VNAME, node->vname, MODEL_KIND_VNAME);
    
    str = lua_tostring(L, -1);
    if (!strcmp(str,KIND_REAL_NAME)) node->kind = REAL_KIND;
#ifdef _HAVE_NUMLUA
    else if (!strcmp(str,KIND_SYMMETRIC_NAME)) {
      node->kind = SYMMETRIC_KIND;
      if (node->Ndims != 2) error (L, "'%s.%s': invalid '%s' wrt. '%s'.",
                                   MODEL_VNAME, node->vname, MODEL_KIND_VNAME,
                                   MODEL_DIM_VNAME);
    }
#endif
    else if (!strcmp(str,KIND_INTEGER_NAME)) node->kind = INTEGER_KIND;
    else error(L, "'%s.%s.%s': Invalid value: '%s'.",
               MODEL_VNAME, node->vname, MODEL_KIND_VNAME, str); 
    lua_pop(L, 1); /*}}}*/
    
    /* .limits */ /*{{{*/
    if (MOREVERBOSE) fprintf(stderr, " --> '%s': %s\n", 
                             node->vname, MODEL_LIMITS_VNAME);
    lua_getfield(L, -1, MODEL_LIMITS_VNAME);
    if (lua_isnil(L,-1)) {
      node->censored = FALSE;
    } else if (node->Ndims != 1 || node->dim != 1) {
      error(L, "'%s.%s' censoring applicable only with univariate"
            " distributions.", MODEL_VNAME, node->vname);
    } else if (!lua_istable(L,-1)) {
      error(L, "'%s.%s.%s' should be a table.", 
            MODEL_VNAME, node->vname, MODEL_LIMITS_VNAME);
    } else {
      /* Everything is fine; continue to set limits */
      node->censored = TRUE;
      /* The minimum value */
      lua_rawgeti(L, -1, 1);
      if (lua_isnil(L,-1))
        node->limits[0] = GRAPHAM_NINF;
      else
        node->limits[0] = luaL_checknumber(L,-1);
      lua_pop(L, 1);
      /* The maximum */
      lua_rawgeti(L, -1, 2);
      if (lua_isnil(L,-1))
        node->limits[1] = GRAPHAM_INF;
      else
        node->limits[1] = luaL_checknumber(L,-1);
      lua_pop(L, 1);
    }
    lua_pop(L, 1);
    /*}}}*/

    /* Allocate memory for value */
    node->value = (double *)malloc(sizeof(double)*node->dim);
    
    /* .init_val */ /*{{{*/
    lua_pushstring(L, MODEL_INIT_VNAME);
    lua_rawget(L, -2);
    if (lua_isnil(L,-1)) {
      if (node->type == CUSTOM_TYPE)
        error(L, "Error: The type of '%s' requires an initial value.",
              node->vname);
      for (j=0;j<node->dim;j++) {
        node->value[j] = 0.0;
      }
    } else {
      if (read_node_value(L, node) != EXIT_SUCCESS)
        error(L, "Cannot read '%s.%s.%s': invalid value.",
              MODEL_VNAME, node->vname, MODEL_INIT_VNAME);
    }
    lua_pop(L, 1); /* init_val */
    if (MOREVERBOSE) {
      fprintf(stderr, " --> '%s': %s [", node->vname, MODEL_INIT_VNAME);
      for (j=0; j<node->dim; j++) 
        fprintf(stderr, "%e,",node->value[j]);
      fprintf(stderr, "]\n");
    }
    /*}}}*/    

    node->in_block = (boolean *)malloc(sizeof(boolean)*node->dim);
    for (j=0; j<node->dim; j++) node->in_block[j] = FALSE;

    if (node->kind == SYMMETRIC_KIND && check_symmetric(node)) 
        error(L, "Error: '%s.%s': incompatible '%s' wrt. '%s'.",
                MODEL_VNAME, node->vname, MODEL_INIT_VNAME, MODEL_KIND_VNAME);

    node->parents = (Node **)malloc(sizeof(Node *)*model->N);
    node->children = (Node **)malloc(sizeof(Node *)*model->N);
    
    lua_pop(L,1); /* value */
    k++;
  } 

  
  if (VERBOSE) 
    fprintf(stderr, " * Going through parents & setting children\n");
  lua_pushnil(L);
  while (lua_next(L, -2) != 0) {
    str = lua_tostring(L, -2);
    k = find_variable(model, str);
    node = &(model->nodes[k]);
    /* .parents */
    lua_getfield(L,-1, MODEL_PARENTS_VNAME);
    if (lua_isnil(L,-1)) {
      node->Nparents = 0;
      node->parents = NULL;
    } else {
      if (!lua_istable(L,-1))
        error(L, "Error: '%s.%s' should be a table.", MODEL_VNAME,
              MODEL_PARENTS_VNAME);
      node->Nparents = lua_objlen(L, -1);
      for (j=0; j<node->Nparents; j++) {
        lua_rawgeti(L, -1, j+1);
        if (!lua_isstring(L,-1)) 
          error(L,"Error: '%s.%s' should be a table of variable names.",
                MODEL_VNAME, MODEL_PARENTS_VNAME);
        
        str = lua_tostring(L, -1);
        i = find_variable(model, str);
        if (i<0) error(L, "Cannot find parent '%s' of variable '%s'\n",
                       str, node->vname);
        for (l=0; l<j; l++) {
          /* Check for duplicates */
          if (node->parents[l] == &(model->nodes[i]))
            error(L, "Duplicate parent '%s' of variable '%s'\n",
                  str, node->vname);
        }
        /* Add parent & children info */
        parent = &(model->nodes[i]);
        node->parents[j] = parent;
        parent->children[parent->Nchildren++] = node;
        
        lua_pop(L,1);
      }
    }
    /* Initialise to zero */
    lua_pop(L, 2); /* value, .parents */
  }
  lua_pop(L, 1); /* model */
  
}

/*}}}*/

/* Try to open a C library */
void* open_clib(const char* fname) { /*{{{*/
  void* handle = NULL;
  char* buf;
  
  /* Try if name is relative to current working dir. */
  buf = (char *)malloc(sizeof(char)*(FILENAME_MAX+PATH_MAX+2));
  if (buf != NULL) {
    if (getcwd(buf, sizeof(char)*PATH_MAX)) {
      strcat(buf, "/");
      strcat(buf, fname);
      handle = dlopen(buf, RTLD_NOW | RTLD_GLOBAL);
    }
  }
  if (handle == NULL) {
    /* Try the default paths */
    handle = dlopen(fname, RTLD_NOW | RTLD_GLOBAL);
  }
  return handle;
} /*}}}*/

/* Read the parameters structure */
void read_parameters(lua_State *L, Model *model, Parameters* para) { /*{{{*/
  /*int k, ind;*/
  char const *fname, *str;
  lua_getglobal(L, PARA_VNAME); 

  /* Simple numeric fields of PARA_VNAME */
  para->niter = read_numeric_field(L, PARA_VNAME, PARA_NITER_VNAME, 
                                   DEFAULT_NITER);
  para->nburn = read_numeric_field(L, PARA_VNAME, PARA_NBURN_VNAME,
                                   DEFAULT_NBURN);
  para->nthin = read_numeric_field(L, PARA_VNAME, PARA_NTHIN_VNAME, 
                                   DEFAULT_NTHIN);
  para->acc_opt1 = read_numeric_field(L, PARA_VNAME, PARA_ACC1_VNAME, 
                                      DEFAULT_ACC_OPT1);
  para->acc_opt2 = read_numeric_field(L, PARA_VNAME, PARA_ACC2_VNAME, 
                                      DEFAULT_ACC_OPT2);
  para->seed = read_numeric_field(L, PARA_VNAME, PARA_SEED_VNAME,
                                  GRAPHAM_NINF);
  para->random_scan = read_numeric_field(L, PARA_VNAME, PARA_RANDOM_SCAN_VNAME, 
                                   DEFAULT_RANDOM_SCAN);
  para->dr_scaling = (read_numeric_field(L, PARA_VNAME, 
                                  PARA_DR_SCALING_VNAME, DEFAULT_DR_SCALING)
                                    == 0.0) ? FALSE : TRUE;
  if (para->dr_scaling<0.0)
    error(L, "Error: %s.%s should be positive.", PARA_VNAME, 
          PARA_DR_SCALING_VNAME);
  
  lua_getglobal(L, PARA_VNAME); 
  lua_getfield(L, -1, PARA_CLIB_VNAME); /*{{{*/
  if (!lua_isnil(L,-1)) {
    if (!lua_isstring(L,-1)) 
    error(L, "'%s.%s' should be a string.", PARA_VNAME, PARA_CLIB_VNAME);
    fname = lua_tostring(L, -1);
    para->clib = open_clib(fname);
    if (para->clib == NULL)
      error(L, "Cannot open dynamic library '%s': %s", fname, dlerror());
  } else {
    para->clib = NULL;
  }
  lua_pop(L,1);
  /*}}}*/

  lua_getfield(L, -1, PARA_ALG_VNAME); /*{{{*/
  if (!lua_isnil(L,-1)) {
    if (!lua_isstring(L,-1)) 
    error(L, "'%s.%s' should be a string.", PARA_VNAME, PARA_ALG_VNAME);
    str = lua_tostring(L, -1);
    if (!strcmp(str,ALG_ASCM_NAME))            para->alg = ALG_ASCM;
    else if (!strcmp(str,ALG_AM_NAME))         para->alg = ALG_AM;
    else if (!strcmp(str,ALG_AMS_NAME))        para->alg = ALG_AMS;
    else if (!strcmp(str,ALG_RBAM_NAME))       para->alg = ALG_RBAM;
    else if (!strcmp(str,ALG_RBAMS_NAME))      para->alg = ALG_RBAMS;
    else if (!strcmp(str,ALG_ASHM_NAME))       para->alg = ALG_ASHM;
    else if (!strcmp(str,ALG_METROPOLIS_NAME)) para->alg = ALG_METROPOLIS;
    else error(L, "Unknown algorithm: '%s'.",str);
  } else {
    para->alg = DEFAULT_ALG;
  }
  lua_pop(L,1);
  /*}}}*/

  lua_getfield(L, -1, PARA_PROPOSAL_VNAME); /*{{{*/
  if (!lua_isnil(L,-1)) {
    if (!lua_isstring(L,-1)) 
    error(L, "'%s.%s' should be a string.", PARA_VNAME, PARA_PROPOSAL_VNAME);
    str = lua_tostring(L, -1);
  } else {
    str = DEFAULT_PROPOSAL;
  }
  if (!strcmp(str,PROPOSAL_GAUSSIAN)) {
    para->rand = &rand_gaussian;
    para->dram_proposal_ratio = &gaussian_ratio;
  } else if (!strcmp(str,PROPOSAL_UNIFORM)) {
    para->rand = &rand_uniform;
    para->dram_proposal_ratio = &uniform_ratio;
  } else if (!strcmp(str,PROPOSAL_STUDENT)) {
    para->rand = &rand_student;
    para->dram_proposal_ratio = &student_ratio;
  } else if (!strcmp(str,PROPOSAL_LAPLACE)) {
    para->rand = &rand_laplace;
    para->dram_proposal_ratio = &laplace_ratio;
  } else {
    error(L, "Unknown proposal distribution: '%s'.",str);
  }
  lua_pop(L,1);
  
  /*}}}*/

  lua_getfield(L, -1, PARA_INIT_VNAME); /*{{{*/
  if (!lua_isnil(L,-1)) {
    if (!lua_isstring(L,-1)) 
    error(L, "'%s.%s' should be a string.", PARA_VNAME, PARA_INIT_VNAME);
    str = lua_tostring(L, -1);
    if (!strcmp(str,INIT_GREEDY_NAME))         para->init = INIT_GREEDY;
    else if (!strcmp(str,INIT_FREEZE_NAME))    para->init = INIT_FREEZE;
    else if (!strcmp(str,INIT_TRAD_NAME))      para->init = INIT_TRAD;
    else error(L, "Unknown initialisation method: '%s'.",str);
  } else {
    para->init = DEFAULT_INIT;
  }
  lua_pop(L,1);
  /*}}}*/

  lua_getfield(L, -1, PARA_BLOCKING_VNAME); /*{{{*/
  if (!lua_isnil(L,-1)) {
    if (!lua_isstring(L,-1)) 
    error(L, "'%s.%s' should be a string.", PARA_VNAME, PARA_BLOCKING_VNAME);
    str = lua_tostring(L, -1);
    if (!strcmp(str,BLOCKING_SC_NAME))        para->blocking = BLOCKING_SC;
    else if (!strcmp(str,BLOCKING_NODE_NAME)) para->blocking = BLOCKING_NODE;
    else if (!strcmp(str,BLOCKING_FULL_NAME)) para->blocking = BLOCKING_FULL;
    else error(L, "Unknown blocking method: '%s'.",str);
  } else {
    para->blocking = DEFAULT_BLOCKING;
  }
  lua_pop(L,1);

/*}}}*/

  lua_getfield(L, -1, PARA_OUTFMT_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isstring(L,-1)) error(L, "'%s.%s' should be a string.",
                                   PARA_VNAME, PARA_OUTFMT_VNAME);
    str = lua_tostring(L,-1);
    if (strcmp(str, OUTFMT_ASCII) == 0) {
      para->output.add_record_outfile = *add_ascii_record;
    } else if (strcmp(str, OUTFMT_BINARY) == 0) {
      para->output.add_record_outfile = *add_binary_record;
    } else {
      error(L, "'%s.%s': invalid value: '%s'.",
            PARA_VNAME, PARA_OUTFMT_VNAME, str);
    }
  } else {
    para->output.add_record_outfile = *add_binary_record;
  }
  lua_pop(L,1);

/*}}}*/

  lua_getfield(L, -1, PARA_OUTF_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isstring(L,-1)) error(L, "'%s.%s' should be a string.",
                                   PARA_VNAME, PARA_OUTF_VNAME);
    fname = lua_tostring(L,-1);
    if (para->output.add_record_outfile == *add_binary_record) {
      para->output.file = open_outfile(L, fname, "wb");
    } else {
      para->output.file = open_outfile(L, fname, "wt");
    }
    para->output.write_file = TRUE;
  } else {
    para->output.write_file = FALSE;
  }
  lua_pop(L,1);
  /*}}}*/

  lua_getfield(L, -1, PARA_OUTADAPTFMT_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isstring(L,-1)) error(L, "'%s.%s' should be a string.",
                                   PARA_VNAME, PARA_OUTADAPTFMT_VNAME);
    str = lua_tostring(L,-1);
    if (strcmp(str, OUTFMT_ASCII) == 0) {
      para->adapt_output.add_record_outfile = *add_ascii_record;
    } else if (strcmp(str, OUTFMT_BINARY) == 0) {
      para->adapt_output.add_record_outfile = *add_binary_record;
    } else {
      error(L, "'%s.%s': invalid value: '%s'.",
            PARA_VNAME, PARA_OUTFMT_VNAME, str);
    }
  } else {
    para->adapt_output.add_record_outfile = *add_binary_record;
  }
  lua_pop(L,1);
  /*}}}*/
  
  lua_getfield(L, -1, PARA_OUTADAPT_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isstring(L,-1)) error(L, "'%s.%s' should be a string.",
                                   PARA_VNAME, PARA_OUTADAPT_VNAME);
    fname = lua_tostring(L,-1);
    if (para->output.add_record_outfile == *add_binary_record) {
      para->adapt_output.file = open_outfile(L, fname, "wb");
    } else {
      para->adapt_output.file = open_outfile(L, fname, "wt");
    }
    para->adapt_output.write_file = TRUE;
  } else {
    para->adapt_output.write_file = FALSE;
  }
  lua_pop(L,1);
  /*}}}*/

  lua_getfield(L, -1, PARA_OUTCFG_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isstring(L,-1)) error(L, "'%s.%s' should be a string.",
                                   PARA_VNAME, PARA_OUTF_VNAME);
    fname = lua_tostring(L,-1);
    para->outcfg = open_outfile(L, fname, "wt");
  } else {
    para->outcfg = NULL;
  }
  lua_pop(L,1);
  
/*}}}*/

  lua_getfield(L, -1, PARA_CLOSE_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isfunction(L, -1)) error(L, "'%s.%s' should be a function.",
                                      PARA_VNAME, PARA_CLOSE_VNAME);
    /* NOTE: pops the function out from the stack! */
    para->close_fun = luaL_ref(L, LUA_REGISTRYINDEX);
    para->is_close_fun = TRUE;
    LUA_CALLS = TRUE;
  } else {
    lua_pop(L, 1);
    para->is_close_fun = FALSE;
  }
  
/*}}}*/

  lua_getfield(L, -1, PARA_CUSTOM_SCALING_FUN_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isfunction(L, -1)) 
      error(L, "Error: '%s.%s' should be a function.",
            PARA_VNAME, PARA_CUSTOM_SCALING_FUN_VNAME);
    /* NOTE: pops the function out from the stack! */
    para->custom_scaling_fun = luaL_ref(L, LUA_REGISTRYINDEX);
    para->is_custom_scaling_fun = TRUE;
  } else {
    lua_pop(L, 1);
    para->is_custom_scaling_fun = FALSE;
  }
  /*}}}*/

  lua_getfield(L, -1, PARA_ADAPT_WEIGHT_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (lua_isfunction(L, -1)) {
      /* NOTE: pops the function out from the stack! */
      para->adapt_weight_fun = luaL_ref(L, LUA_REGISTRYINDEX);
      para->is_adapt_weight_fun = TRUE;
    } else if (lua_isnumber(L, -1)) {
      para->is_adapt_weight_fun = FALSE;
      para->adapt_weight_exp = lua_tonumber(L, -1);
      if (para->adapt_weight_exp < 0.0) 
        error(L, "Error: '%s.%s' should be non-negative.",
              PARA_VNAME, PARA_ADAPT_WEIGHT_VNAME);
      lua_pop(L, 1);
    } else {
        error(L, "Error: '%s.%s' should be a number or a function.",
              PARA_VNAME, PARA_ADAPT_WEIGHT_VNAME);
    }
  } else {
    lua_pop(L, 1);
    para->is_adapt_weight_fun = FALSE;
    para->adapt_weight_exp = DEFAULT_ADAPT_WEIGHT;
  }
  /*}}}*/

  lua_getfield(L, -1, PARA_ADAPT_WEIGHT_SC_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (lua_isfunction(L, -1)) {
      /* NOTE: pops the function out from the stack! */
      para->adapt_weight_sc_fun = luaL_ref(L, LUA_REGISTRYINDEX);
      para->is_adapt_weight_sc_fun = TRUE;
    } else if (lua_isnumber(L, -1)) {
      para->is_adapt_weight_sc_fun = FALSE;
      para->adapt_weight_sc_exp = lua_tonumber(L, -1);
      if (para->adapt_weight_sc_exp < 0.0) 
        error(L, "Error: '%s.%s' should be non-negative.",
              PARA_VNAME, PARA_ADAPT_WEIGHT_SC_VNAME);
      lua_pop(L, 1);
    } else {
        error(L, "Error: '%s.%s' should be a number or a function.",
              PARA_VNAME, PARA_ADAPT_WEIGHT_SC_VNAME);
    }
  } else {
    lua_pop(L, 1);
    para->is_adapt_weight_sc_fun = FALSE;
    para->adapt_weight_sc_exp = DEFAULT_ADAPT_WEIGHT_SC;
  }
  /*}}}*/

  lua_getfield(L, -1, PARA_MIX_WEIGHT_VNAME); /*{{{*/
  if (!lua_isnil(L, -1)) {
    if (!lua_isfunction(L, -1)) 
      error(L, "Error: '%s.%s' should be a function.",
            PARA_VNAME, PARA_MIX_WEIGHT_VNAME);
    /* NOTE: pops the function out from the stack! */
    para->mix_weight_fun = luaL_ref(L, LUA_REGISTRYINDEX);
    para->is_mix_weight_fun = TRUE;
  } else {
    lua_pop(L, 1);
    para->is_mix_weight_fun = FALSE;
  }
  /*}}}*/
  
  lua_pop(L,1); /* PARA_NAME */

}

/*}}}*/

/* Read the functional */
void read_functional(lua_State *L, Model* model, Parameters* para, 
                     Functional* functional) { /*{{{*/
  int j, k, N;
  const char* fname;
  Node* node;
  
  lua_getglobal(L, FUNC_VNAME);
  if (lua_isnil(L, -1)) { /*{{{*/
    functional->evaluate = FALSE;
    functional->called = 0;
    lua_pop(L, 1); /* nil */
    /*}}}*/
  } else if (lua_isfunction(L, -1)) { /*{{{*/
    functional->evaluate = TRUE;
    /* this pops the stack! */
    functional->fun = luaL_ref(L, LUA_REGISTRYINDEX);
    functional->called = 0;
    functional->is_cfun = FALSE;
    LUA_CALLS = TRUE;
    /*}}}*/
  } else if (lua_istable(L, -1)) { /*{{{*/
    /* a c function */
    functional->called = 0;
    functional->is_cfun = TRUE;
    functional->evaluate = TRUE;
    functional->mean_N = read_numeric_field(L, FUNC_VNAME, 
                                            FUNC_DIM_VNAME, -1);
    if (functional->mean_N < 0) 
      error(L,"'%s.%s' should be a positive integer.", 
            FUNC_VNAME, FUNC_DIM_VNAME);
    
    lua_getfield(L, -1, FUNC_NAME_VNAME); 
    if (!lua_isstring(L,-1)) 
      error(L, "'%s.%s' should be a string.", FUNC_VNAME, FUNC_NAME_VNAME);
    fname = lua_tostring(L, -1);
    
    functional->cfun = NULL;
    /* Match from user-defined lib */
    if (para->clib != NULL) {
      *(void **)(&functional->cfun) = dlsym(para->clib, fname);
    }
    /* Match from grapham_flib */
    if (functional->cfun == NULL) {
      for (j=0; ; j++) {
        if (grapham_flib[j].name == NULL) {
          break;
        } else {
          if (strcmp(grapham_flib[j].name, fname) == 0) {
            functional->cfun = grapham_flib[j].c_fun;
          }
        }
      }
    } 
    if (functional->cfun == NULL)
      error(L, "'%s.%s': Cannot find matching function.",
            FUNC_VNAME, FUNC_NAME_VNAME);
    lua_pop(L, 1); /* name */
    
    lua_getfield(L, -1, FUNC_ARGS_VNAME);
    if (!lua_istable(L, -1))
      error(L, "'%s.%s' should be a table.",FUNC_VNAME, FUNC_ARGS_VNAME);
    
    N = lua_objlen(L, -1);
    functional->cfun_Nargs = N;
    functional->cfun_args = (const double **)malloc(sizeof(double *)*N);
    functional->cfun_lengths = (int *)malloc(sizeof(int)*N);
    for (j=0; j<N; j++) {
      lua_rawgeti(L, -1, j+1);
      fname = luaL_checkstring(L, -1);
      k = find_variable(model, fname);
      if (k<0)
        error(L, "Cannot find variable '%s' in model.",
              fname);
      node = &(model->nodes[k]);
      functional->cfun_args[j] = node->value;
      functional->cfun_lengths[j] = node->dim;
      lua_pop(L, 1); /* field */
    }
    lua_pop(L, 1); /* args */
    
    /*}}}*/
  } else {
    error(L,"'%s' must be a function or a table.", FUNC_VNAME);
  }
} /*}}}*/


/* Write the adaptation parameters (Cholesky factor & scaling) for each block,
 * as well as the block structure, to an output file. */
void write_outcfg(FILE* ofile, Model* model) { /*{{{*/
  int i,j,k, ind;
  Block* block;
  Node* node;
  
  fprintf(ofile, "%s.%s = {\n", PARA_VNAME, PARA_BLOCKS_VNAME);
  for (k=0; k<model->Nblocks; k++) {
    block = model->blocks[k];
    fprintf(ofile, "  {");
    for (i=0; i<block->dim; i++) {
      node = block->components[i].node;
      ind = block->components[i].index;
      if (node->type == NUMBER_TYPE) {
        fprintf(ofile, "\"%s\",", node->vname);
      } else {
        fprintf(ofile, "\"%s[%i]\",", node->vname, ind+1);
      }
    }
    fprintf(ofile, "},\n");
  }
  fprintf(ofile, "}\n\n");

  fprintf(ofile, "%s.%s = {\n", PARA_VNAME, PARA_BLOCKS_CHOL_VNAME);
  for (k=0; k<model->Nblocks; k++) {
    block = model->blocks[k];
    fprintf(ofile, "  {\n");
    for (i=0; i<block->dim; i++) {
      fprintf(ofile, "    {");
      for (j=0; j<block->dim; j++) {
        fprintf(ofile, OUTCFG_NUMFORMAT, block->adapt_chol[IND(i,j,block->dim)]);
      }
      fprintf(ofile, "},\n");
    }
    fprintf(ofile, "  },\n");
  }
  fprintf(ofile, "}\n\n");
  
  fprintf(ofile, "%s.%s = {\n", PARA_VNAME, PARA_BLOCKS_SC_VNAME);
  for (k=0; k<model->Nblocks; k++) {
    block = model->blocks[k];
    fprintf(ofile, "  ");
    fprintf(ofile, OUTCFG_NUMFORMAT, block->scaling);
    fprintf(ofile, "\n");
  }
  fprintf(ofile, "}\n\n");
  fclose(ofile);  
} /*}}}*/

/* Check that the graph is acyclic */
int mark_connected_children(Node* node, Node* root) { /*{{{*/
  boolean cn;
  int k;
  if (!node->connected) {
    node->connected = TRUE;
    for (k=0; k<node->Nchildren; k++) {
      if (node->children[k] == root) {
        return EXIT_FAILURE;
      } else {
        cn = mark_connected_children(node->children[k], root);
        if (cn != EXIT_SUCCESS) return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
} /*}}}*/
int check_dag(Model* model) { /*{{{*/
  int k, j;
  boolean connected = FALSE;
  for (j=0; j<model->N; j++) {
    /* This may be far from optimal, but does not consume memory... */
    for (k=0; k<model->N; k++) {
      model->nodes[k].connected = FALSE;
    }
    connected = mark_connected_children(&(model->nodes[j]), &(model->nodes[j]));
    if (connected) break;
  }
  return connected;
}

/*}}}*/

/* Check that the graph is connected */
void mark_connected_neighbours(Node* node) { /*{{{*/
  int k;
  if (!node->connected) {
    node->connected = TRUE;
    for (k=0; k<node->Nchildren; k++) {
      mark_connected_neighbours(node->children[k]);
    }
    for (k=0; k<node->Nparents; k++) {
      mark_connected_neighbours(node->parents[k]);
    }
  }
  return;
} /*}}}*/
void check_connectivity(Model* model) { /*{{{*/
  int k, root_ind = 0;
  boolean not_connected = FALSE;
  for (k=0; k<model->N; k++) {
    model->nodes[k].connected = FALSE;
    /* Make sure the root node not instantiated, if any! */
    if (!model->nodes[k].instantiated) root_ind = k;
  }
  mark_connected_neighbours(&(model->nodes[root_ind]));
  for (k=0; k<model->N; k++) {
    if (!model->nodes[k].connected
        && !model->nodes[k].instantiated) not_connected = TRUE;
  }
  if (not_connected) {
    fprintf(stderr, "Warning: the following nodes:\n  ");
    for (k=0; k<model->N; k++) {
      if (!model->nodes[k].connected) 
      fprintf(stderr, " '%s'", model->nodes[k].vname);
    }
    fprintf(stderr, "\nare not connected with:\n  ");
    for (k=0; k<model->N; k++) {
      if (model->nodes[k].connected) 
      fprintf(stderr, " '%s'", model->nodes[k].vname);
    }
    fprintf(stderr, "\n\n");
  }
} /*}}}*/

/* Just go through the predefined blocks, check sanity
 * and mark each variable 'in_block' flag, if it is
 * already in a block. */
void mark_blocks_used(lua_State *L, Model* model) { /*{{{*/
  Node* node;
  int k, j, i, ind, ind_e, Nblocks;
  const char* str;
  
  if (MOREVERBOSE) fprintf(stderr, " * Going through sampling blocks.\n");
  lua_getglobal(L, PARA_VNAME);
  lua_getfield(L, -1, PARA_BLOCKS_VNAME);
  if (lua_isnil(L, -1)) {
    if (MOREVERBOSE) fprintf(stderr, " --> Found no blocks.\n");
    lua_pop(L, 1); /* nil */
    lua_newtable(L);
    lua_setfield(L, -2, PARA_BLOCKS_VNAME);
    lua_getfield(L, -1, PARA_BLOCKS_VNAME);
    Nblocks = 0;
  } else {
    if (!lua_istable(L, -1)) {
      error(L, "'%s.%s' should be a table of tables of variable names.", 
            PARA_VNAME, PARA_BLOCKS_VNAME);
    }
    /* Go through the existing blocks & check for validity & 
     * just mark what nodes are used */
    Nblocks = lua_objlen(L,-1);
    if (MOREVERBOSE) fprintf(stderr, " --> Found %i blocks.\n", Nblocks);
    for (k=0; k<Nblocks; k++) {
      lua_rawgeti(L, -1, k+1);
      if (!lua_istable(L, -1)) {
        error(L, "'%s.%s' should be a table of tables of variable names.", 
              PARA_VNAME, PARA_BLOCKS_VNAME);
      }
      for (j=0; j<lua_objlen(L,-1); j++) {
        lua_rawgeti(L, -1, j+1);
        if (!lua_isstring(L, -1)) {
          error(L, "'%s.%s' should be a table of tables of variable names.", 
                PARA_VNAME, PARA_BLOCKS_VNAME);
        }
        str = lua_tostring(L, -1);
        ind = find_variable_index(model, str, &ind_e);
        if (ind<0) error(L, "'%s.%s': Unknown variable '%s'.", 
                         PARA_VNAME, PARA_BLOCKS_VNAME, str);
        node = &(model->nodes[ind]);
        if (node->instantiated)
        error(L, "'%s.%s': Variable '%s' in instantiated.", 
              PARA_VNAME, PARA_BLOCKS_VNAME, str);
        if (ind_e<0) {
          for (i=0; i<node->dim; i++) {
            if (node->in_block[i])
            error(L, "'%s.%s': Variable '%s[%i]' in two blocks.", 
                  PARA_VNAME, PARA_BLOCKS_VNAME, str, i);
            node->in_block[i] = TRUE;
          }
        } else {
          if (node->in_block[ind_e])
          error(L, "'%s.%s': Variable '%s[%i]' in two blocks.", 
                PARA_VNAME, PARA_BLOCKS_VNAME, str, ind_e);
          node->in_block[ind_e] = TRUE;
        }
        lua_pop(L, 1); /* str */
      }
      lua_pop(L, 1); /* j:th block (table) */
    } 
  }
  lua_pop(L, 2); /* para.blocks */
} /*}}}*/

/* Add a block for each uninstantiated variable that is
 * not currently attached to any block */
void fill_blocks(lua_State *L, Model* model) { /*{{{*/
  Node* node;
  int k, j, n, Nblocks;
  char* tmpstr;
  
  mark_blocks_used(L, model);
  
  lua_getglobal(L, PARA_VNAME);
  lua_getfield(L, -1, PARA_BLOCKS_VNAME);
  Nblocks = lua_objlen(L, 1);
  
  if (MOREVERBOSE) 
  fprintf(stderr, " * Filling in the rest variables as individual blocks.\n");
  /* The actual fill */
  for (k=0; k<model->N; k++) {
    node = &(model->nodes[k]);
    if (node->instantiated) continue;
    n = 0;
    lua_newtable(L);
    for (j=0; j<node->dim; j++) {
      if (node->in_block[j]) continue;
      tmpstr = (char *)malloc(sizeof(char)*strlen(node->vname)+30);
      sprintf(tmpstr, "%s[%i]", node->vname, j+1);
      lua_pushstring(L, tmpstr);
      lua_rawseti(L, -2, ++n);
      free(tmpstr);
    }
    if (n==0) {
      lua_pop(L, 1); /* discard the empty table */
    } else {
      lua_rawseti(L, -2, ++Nblocks); /* newtable->para.blocks */
    }
  }
  lua_pop(L, 2); /* para, .blocks */
}
/*}}}*/

/* Add all the rest free variables to a single block */
void fill_single_block(lua_State *L, Model* model) { /*{{{*/
  Node* node;
  int k, j, n, Nblocks;
  char* tmpstr;
  
  mark_blocks_used(L, model);
  
  lua_getglobal(L, PARA_VNAME);
  lua_getfield(L, -1, PARA_BLOCKS_VNAME);
  Nblocks = lua_objlen(L, 1);
  
  if (MOREVERBOSE) 
  fprintf(stderr, " * Filling in the rest variables as one block.\n");
  /* The actual fill */
  n = 0;
  lua_newtable(L);
  for (k=0; k<model->N; k++) {
    node = &(model->nodes[k]);
    if (node->instantiated) continue;
    for (j=0; j<node->dim; j++) {
      if (node->in_block[j]) continue;
      tmpstr = (char *)malloc(sizeof(char)*strlen(node->vname)+30);
      sprintf(tmpstr, "%s[%i]", node->vname, j+1);
      lua_pushstring(L, tmpstr);
      lua_rawseti(L, -2, ++n);
      free(tmpstr);
    }
  }
  if (n==0) {
    lua_pop(L, 1); /* discard the empty table */
  } else {
    lua_rawseti(L, -2, ++Nblocks); /* newtable->para.blocks */
  }
  lua_pop(L, 2); /* para, .blocks */
}
/*}}}*/

/* Add each component not instantiated or already in a sampling block
 * to a separate block */
void fill_single_components(lua_State *L, Model* model) { /*{{{*/
  Node* node;
  int k, j, Nblocks;
  char* tmpstr;
  
  mark_blocks_used(L, model);
  
  lua_getglobal(L, PARA_VNAME);
  lua_getfield(L, -1, PARA_BLOCKS_VNAME);
  Nblocks = lua_objlen(L, 1);
  
  if (MOREVERBOSE) 
  fprintf(stderr, " * Filling in the rest components as individual blocks.\n");
  /* The actual fill */
  for (k=0; k<model->N; k++) {
    node = &(model->nodes[k]);
    if (node->instantiated) continue;
    for (j=0; j<node->dim; j++) {
      if (node->in_block[j]) continue;
      tmpstr = (char *)malloc(sizeof(char)*strlen(node->vname)+30);
      lua_newtable(L);
      sprintf(tmpstr, "%s[%i]", node->vname, j+1);
      lua_pushstring(L, tmpstr);
      lua_rawseti(L, -2, 1);
      free(tmpstr);
      lua_rawseti(L, -2, ++Nblocks);
    }
  }
  lua_pop(L, 2); /* para, .blocks */
}
/*}}}*/

/* Set up the sampling blocks, according to para.block. */
void read_blocks(lua_State *L, Model* model) { /*{{{*/
  int j, i, k, m, ind, ind_e, Nblocks, Nvars, max_dim = 0;
  double num = 0.0, *work;
  const char* str;
  Block* block; Node* node;
  
  lua_getglobal(L, PARA_VNAME);
  lua_getfield(L, -1, PARA_BLOCKS_VNAME);
  
  if (MOREVERBOSE)
  fprintf(stderr, " --> Reading '%s.%s'.\n", 
          PARA_VNAME, PARA_BLOCKS_VNAME);
  
  Nblocks = lua_objlen(L, -1);
  model->blocks = (Block **)malloc(sizeof(Block*)*Nblocks);
  for (k=0; k<Nblocks; k++) { 
    block = (Block *)malloc(sizeof(Block));
    /* Through blocks */
    lua_rawgeti(L, -1, k+1);
    if (!lua_istable(L, -1)) {
      error(L, "'%s.%s' should be a table of tables of variable names.", 
            PARA_VNAME, PARA_BLOCKS_VNAME);
    }
    block->dim = 0;
    block->Nchildren = 0;
    block->variables = (Node **)malloc(sizeof(Node *)*model->N);
    block->components = (Component *)malloc(sizeof(Component)*model->dim);
    block->children = (Node **)malloc(sizeof(Node *)*model->N);
    
    /* Determine the components in the block */
    for (j=0; j<lua_objlen(L, -1); j++) {
      lua_rawgeti(L, -1, j+1);
      str = lua_tostring(L, -1);
      ind = find_variable_index(model, str, &ind_e);
      lua_pop(L, 1);
      node = &(model->nodes[ind]);
      if (ind < 0)
        error(L, "'%s.%s': Cannot find '%s'", 
            PARA_VNAME, PARA_BLOCKS_VNAME, str);
      if (ind_e < 0) {
        for (i=0; i<node->dim; i++) {
          block->components[block->dim].index = i;
          block->components[block->dim].node = node;
          block->dim++;
        }
      } else {
        block->components[block->dim].index = ind_e;
        block->components[block->dim].node = node;
        block->dim++;
      }
    }
    block->components = (Component *)
    realloc(block->components, sizeof(Component)*block->dim);
    
    /* The default scaling */
    block->init_scaling = DEFAULT_SCALING(block->dim);
    
    max_dim = MAX(block->dim, max_dim);
    
    /* The accumulated mean and Cholesky factor */
    block->adapt_mean = (double *)malloc(sizeof(double)*block->dim);
    block->adapt_chol = (double *)malloc(sizeof(double)*(block->dim)*(block->dim));
    block->init_chol  = (double *)malloc(sizeof(double)*(block->dim)*(block->dim));

    /* Initialise counters of acceptance and rejection */
    block->accepted = 0;
    block->rejected = 0;
    
    /* Set the default value for Cholesky factor */
    set_identity_matrix(block->init_chol, block->dim);

    /* Initialise mean to current value */
    for (j=0; j<block->dim; j++) {
      node = block->components[j].node;
      ind = block->components[j].index;
      block->adapt_mean[j] = node->value[ind];
    }
    
    /* Set variables & children */
    Nvars = 0;
    for (j=0; j<block->dim; j++) {
      node = block->components[j].node;
      for (i=0; i<j; i++) {
        if (block->components[i].node == node) goto next_component;
      }
      block->variables[Nvars++] = node;
      for (i=0; i<node->Nchildren; i++) {
        for (m=0; m<block->Nchildren; m++) {
          if (node->children[i] == block->children[m]) goto next_child;
        }
        block->children[block->Nchildren++] = node->children[i];
        next_child:
        ;
      }
      next_component:
      ;
    }
    lua_pop(L, 1);
    block->Nvars = Nvars;
    block->children = (Node **)realloc(block->children, 
                                       sizeof(Node *)*block->Nchildren);
    block->variables = (Node **)realloc(block->variables, 
                                        sizeof(Node *)*block->Nvars);
    
    model->blocks[k] = block;
  }

  /* Final touch */
  model->Nblocks = Nblocks;
  
  lua_pop(L, 1); /* .blocks */
  
  /* The temporary storage: */
  work = (double *)malloc(sizeof(double)*6*max_dim);
  for (k=0; k<Nblocks; k++) {
    block = model->blocks[k];
    /*block->value = (double *)malloc(sizeof(double)*block->dim);
    block->old_value = (double *)malloc(sizeof(double)*block->dim);
    block->tmp_value = (double *)malloc(sizeof(double)*(block->dim)*3);*/
    block->proposal_value = work;
    block->old_value = &(work[block->dim]);
    block->rand_value = &(work[2*block->dim]);
    block->tmp_value = &(work[3*block->dim]);
  }
  
  if (MOREVERBOSE)
  fprintf(stderr, " --> Reading '%s.%s'.\n", 
          PARA_VNAME, PARA_BLOCKS_CHOL_VNAME);
  /* Load predefined Cholesky factors */
  lua_getfield(L, -1, PARA_BLOCKS_CHOL_VNAME);
  if (!lua_isnil(L, -1)) {
    if (!lua_istable(L, -1)) 
    error(L, "Error: '%s.%s' should be a table.", PARA_VNAME, 
          PARA_BLOCKS_CHOL_VNAME);
    for (k=0; k<Nblocks; k++) {
      block = model->blocks[k];
      lua_rawgeti(L, -1, k+1);
      if (lua_isnil(L, -1)) goto next_chol;
      if (!lua_istable(L, -1))
        error(L, "Error: '%s.%s[%i]' should be a table.", PARA_VNAME, 
              PARA_BLOCKS_CHOL_VNAME, k+1);
      for (j=0; j<block->dim; j++)  {
        lua_rawgeti(L, -1, j+1);
        if (!lua_istable(L, -1))
          error(L, "Error: '%s.%s[%i][%i]' should be a table.", PARA_VNAME, 
                PARA_BLOCKS_CHOL_VNAME, k+1, j+1);
        for (i=0; i<block->dim; i++) {
          lua_rawgeti(L, -1, i+1);
          if (lua_isnil(L, -1)) {
            num = 0.0;
          } else if (lua_isnumber(L, -1)) {
            num = lua_tonumber(L, -1);
          } else {
            error(L, "Error: '%s.%s[%i][%i][%i]' should be a number.", 
                  PARA_VNAME, PARA_BLOCKS_CHOL_VNAME, k+1, j+1, i+1);
          }
          block->init_chol[IND(j,i,block->dim)] = num;
          lua_pop(L, 1);
        }
        lua_pop(L, 1);
      }
      next_chol:
      lua_pop(L, 1); /* (k+1)'th element */
    }
  }
  lua_pop(L, 1); /* .block_chol */

  if (MOREVERBOSE)
  fprintf(stderr, " --> Reading '%s.%s'.\n", 
          PARA_VNAME, PARA_BLOCKS_CHOL_SC_VNAME);
  /* Load and apply additional Cholesky scaling factors */
  lua_getfield(L, -1, PARA_BLOCKS_CHOL_SC_VNAME);
  if (!lua_isnil(L, -1)) {
    if (!lua_istable(L, -1) && !lua_isnumber(L, -1)) 
    error(L, "Error: '%s.%s' should be a table or number.", PARA_VNAME, 
          PARA_BLOCKS_CHOL_SC_VNAME);
    for (k=0; k<Nblocks; k++) {
      block = model->blocks[k];
      if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, k+1);
        if (!lua_isnil(L, -1)) {
          if (!lua_isnumber(L, -1))
          error(L, "Error: '%s.%s[%i]' should be a number.", PARA_VNAME, 
                PARA_BLOCKS_CHOL_SC_VNAME, k+1);
          num = lua_tonumber(L, -1);
        } else {
          num = 1.0;
        }
        lua_pop(L, 1); /* (k+1)'th element */
      } else {
        num = lua_tonumber(L, -1);
      }
      scale_triu(block->init_chol, num, block->dim);
    }
  }
  lua_pop(L, 1); /* .blocks_chol_scaling */
  
  if (MOREVERBOSE)
  fprintf(stderr, " --> Reading '%s.%s'.\n", PARA_VNAME, PARA_BLOCKS_SC_VNAME);
  /* Load predefined scaling factors */
  lua_getfield(L, -1, PARA_BLOCKS_SC_VNAME);
  if (!lua_isnil(L, -1)) {
    if (!lua_istable(L, -1) && !lua_isnumber(L, -1)) 
    error(L, "Error: '%s.%s' should be a table or number.", PARA_VNAME, 
          PARA_BLOCKS_SC_VNAME);
    for (k=0; k<Nblocks; k++) {
      block = model->blocks[k];
      if (lua_istable(L, -1)) {
        lua_rawgeti(L, -1, k+1);
        if (!lua_isnil(L, -1)) {
          if (!lua_isnumber(L, -1))
          error(L, "Error: '%s.%s[%i]' should be a number.", PARA_VNAME, 
                PARA_BLOCKS_SC_VNAME, k+1);
          block->init_scaling = lua_tonumber(L, -1);
        }
        lua_pop(L, 1); /* (k+1)'th element */
      } else {
        block->init_scaling = lua_tonumber(L, -1);
      }
    }
  }
  lua_pop(L, 1); /* .blocks_sc */
  lua_pop(L, 1); /* para */
  
  
  for (k=0; k<Nblocks; k++) {
    block = model->blocks[k];
    memcpy(block->adapt_chol, block->init_chol, 
           sizeof(double)*(block->dim)*(block->dim));
    block->scaling = block->init_scaling;
  }
}
/*}}}*/

