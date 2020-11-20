/* -*- mode: C; mode: fold -*- */
/*
 * src/lua_tools.c
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

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include <string.h>

#ifdef _HAVE_NUMLUA
#include <luamatrix.h>
#endif

#include "grapham_math.h"
#include "grapham_types.h"

#include "lua_tools.h"

/* General error function */
void error (lua_State *L, const char *fmt, ...) { /*{{{*/
  va_list argp;
  va_start(argp, fmt);
  vfprintf(stderr, fmt, argp);
  fprintf(stderr, "\n");
  va_end(argp);
  lua_close(L);
  exit(EXIT_FAILURE);
}

/*}}}*/

Var_type lua_gettype(lua_State* L, int ind, int* dim, int size[2], 
                     int* Ndims) { /*{{{*/
#ifdef _HAVE_NUMLUA
  int match;
  lua_Matrix *M;
#endif
  Var_type type = UNKNOWN_TYPE;
  if (lua_isnumber(L, ind)) {
    dim[0] = size[0] = Ndims[0] = 1;
    type = NUMBER_TYPE;
  } else if (lua_istable(L, ind)) {
    size[0] = lua_objlen(L, ind);
    if (size[0] > 0) {
      lua_rawgeti(L, ind, 1);
      if (lua_istable(L, -1)) {
        Ndims[0] = 2;
        size[1] = lua_objlen(L, -1);
        dim[0] = size[0]*size[1];
      } else {
        Ndims[0] = 1;
        dim[0] = size[0];
      }
      lua_pop(L, 1);
      type = VECTOR_TYPE;
    }
  } else if (lua_isuserdata(L, ind)) {
#ifdef _HAVE_NUMLUA
    lua_getmetatable(L, ind);
    lua_getfield(L, LUA_REGISTRYINDEX, LUAMATRIX_MT);
    match = lua_rawequal(L, -1, -2);
    lua_pop(L, 2);
    if (match) {
    /* This does not work, if not called from Lua! */
    /*if (matrix_istype(L, ind)) {*/
      M = (lua_Matrix *)lua_touserdata(L, ind);
      Ndims[0] = M->dim;
      type = MATRIX_TYPE;
      if (M->dim == 1) {
        dim[0] = size[0] = M->size;
      } else if (M->dim == 2) {
        size[0] = M->size;
        size[1] = M->level[0]->size;
        dim[0] = size[0]*size[1];
      } else {
        type = UNKNOWN_TYPE;
      }
    } else 
#endif
    {
    }
  }
  return type;
} /*}}}*/

int lua_getlength(lua_State* L, int ind) { /*{{{*/
  /* Try to resolve length of a Lua object at position 'ind'.
   * The options are:
   * - A table, then return length
   * - Else, try to call a field "size"
   */
  int N = -1;
#ifdef _HAVE_NUMLUA
  lua_Matrix *M;
#endif
  if (lua_istable(L, ind)) {
    N = lua_objlen(L, ind);
#ifdef _HAVE_NUMLUA
  } else if (matrix_istype(L, ind)) {
    M = (lua_Matrix*)lua_touserdata(L, ind);
    N = M->size;
#endif
  } else {
    lua_getfield(L, ind, "size");
    if (lua_isfunction(L, -1)) {
      lua_pushvalue(L, ind);
      lua_call(L, 1, 1);
      N = luaL_checknumber(L, -1);
      lua_pop(L, 1); /* the returned value */
    }
  }
  return N;
}

/*}}}*/

int lua_getvector(lua_State* L, int ind, double* x, int Nx) { /*{{{*/
#ifdef _HAVE_NUMLUA
  lua_Matrix *M;
#endif
  int k;
  if (lua_istable(L, ind)) {
    for (k=0; k<Nx; k++) {
      lua_rawgeti(L, ind, k+1); x[k] = lua_tonumber(L, -1); lua_pop(L, 1);
    }
#ifdef _HAVE_NUMLUA
  } else if (matrix_istype(L, ind)) {
    M = (lua_Matrix *)lua_touserdata(L, ind);
    memcpy(x, M->data, sizeof(double)*Nx);
    /*for (k=0; k<Nx; k++) printf("%e,",data[k]); printf("\n");*/
#endif
  } else {
    for (k=0; k<Nx; k++) {
      lua_pushnumber(L, k+1); lua_gettable(L, ind);
      x[k] = lua_tonumber(L, -1); lua_pop(L, 1);
    }
  }
  return 1;
} /*}}}*/

/* Read lower-triangular array from Lua to upper-triangular
 * aray in C. */
int lua_gettril(lua_State* L, int ind, double* x, int Nx) { /*{{{*/
#ifdef _HAVE_NUMLUA
  lua_Matrix *M;
#endif
  int k, j;
  if (lua_istable(L, ind)) {
    for (k=0; k<Nx; k++) {
      lua_rawgeti(L, ind, k+1); 
      for (j=0; j<=k; j++) {
        lua_rawgeti(L, -1, j+1); 
        x[IND(j,k,Nx)] = lua_tonumber(L, -1);
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
#ifdef _HAVE_NUMLUA
  } else if (matrix_istype(L, ind)) {
    M = (lua_Matrix*)lua_touserdata(L, ind);
    luaL_argcheck(L, 
                  (M->dim == 2) && (M->size == Nx) 
                  && (M->level[0]->size == Nx),
                  ind, "square matrix expexted.");
    if (istype_upper(M)) {
      for (k=0; k<Nx; k++) 
        for (j=0; j<=k; j++) 
          x[IND(j,k,Nx)] = M->data[IND(j,k,Nx)];
    } else {
      for (k=0; k<Nx; k++) 
        for (j=0; j<=k; j++) 
          x[IND(j,k,Nx)] = M->data[IND(k,j,Nx)];
    }
#endif
  } else {
    for (k=0; k<Nx; k++) {
      lua_pushnumber(L, k+1); lua_gettable(L, ind);
      for (j=0; j<=k; j++) {
        lua_pushnumber(L, j+1); lua_gettable(L, -2);
        x[IND(j,k,Nx)] = lua_tonumber(L, -1); 
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
  }
  return 1;
} /*}}}*/
