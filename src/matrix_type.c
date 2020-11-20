/* -*- mode: C; mode: fold -*- */
/*
 * src/matrix_type.c
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

#include "matrix_type.h"
#include <stdlib.h>
#include <string.h>
#include <luamatrix.h>

void matrix_set_single_value(lua_State* L, double value, 
                             Node* node, const int ind) { /*{{{*/
  int sz;
  lua_Matrix* M;
  
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  M = (lua_Matrix *)lua_touserdata(L, -1);
  if (node->Ndims == 1) {
    M->data[ind] = value;
  } else { /* Ndims == 2 */
    M->data[ind] = value;
    if (node->kind == SYMMETRIC_KIND) {
      sz = node->size[0];
      M->data[(ind%sz)*sz + (ind/sz)] = value;
    }
  }
  lua_pop(L,1);
} /*}}}*/

void matrix_set_node(lua_State* L, Node* node) { /*{{{*/
  lua_Matrix* M;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  M = (lua_Matrix *)lua_touserdata(L, -1);
  memcpy(M->data, node->value, sizeof(double)*node->dim);
  lua_pop(L, 1);
}
/*}}}*/

void matrix_set_node_block(lua_State* L, Node* node, Block* block) { /*{{{*/
  int sz, k, ind;
  lua_Matrix* M;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  M = (lua_Matrix *)lua_touserdata(L, -1);
  for (k=0; k<block->dim; k++) {
    if (block->components[k].node != node) continue;
    ind = block->components[k].index;
    M->data[ind] = block->value[k];
    if (node->kind == SYMMETRIC_KIND) {
      sz = node->size[0];
      M->data[(ind%sz)*sz + (ind/sz)] = block->value[k];
    }
  }
  lua_pop(L, 1); /* node->lua_value */
} /*}}}*/

int matrix_read_node(lua_State* L, Node* node) { /*{{{*/
  lua_Matrix* M = (lua_Matrix *)lua_touserdata(L, -1);
  if (M == NULL) return EXIT_FAILURE;
  if (node->Ndims == 1) {
    if (M->dim != 1 || M->size < node->dim) return -1;
    memcpy(node->value, M->data, sizeof(double)*node->dim);
  } else {
    if (M->dim != 2 || M->size*M->level[0]->size < node->dim) return -1;
    memcpy(node->value, M->data, sizeof(double)*node->dim);
  }
  return EXIT_SUCCESS;
} /*}}}*/

void matrix_create(lua_State* L, Node* node) { /*{{{*/
  lua_Matrix* M;
  /* Create the matrix */
  lua_getglobal(L, "matrix");
  if (node->Ndims == 1) {
    lua_pushnumber(L, node->dim);
    lua_call(L, 1, 1);
  } else {
    lua_pushnumber(L, node->size[0]);
    lua_pushnumber(L, node->size[1]);
    lua_call(L, 2, 1);
  }
  /* Set initial values */
  M = (lua_Matrix *)lua_touserdata(L, -1);
  memcpy(M->data, node->value, sizeof(double)*node->dim);
  lua_setglobal(L, node->vname);
} /*}}}*/
