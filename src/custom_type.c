/* -*- mode: C; mode: fold -*- */
/*
 * src/custom_type.c
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
#include <stdlib.h>
#include "custom_type.h"

void custom_set_single_value(lua_State* L, double value, 
                             Node* node, const int ind) { /*{{{*/
  int i, j, sz;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  if (node->Ndims == 1) {
    lua_pushnumber(L, ind+1);
    lua_pushnumber(L, value);
    lua_settable(L, -3);
  } else {
    sz = node->size[0];
    i = (ind%sz) + 1;
    j = (ind/sz) + 1;
    lua_pushnumber(L, i);
    lua_gettable(L, -2);
    lua_pushnumber(L, j);
    lua_pushnumber(L, value);
    lua_settable(L, -3);
    lua_pop(L, 1);
    if (node->kind == SYMMETRIC_KIND && i != j) {
      lua_pushnumber(L, j + 1);
      lua_gettable(L, -2);
      lua_pushnumber(L, i + 1);
      lua_pushnumber(L, value);
      lua_settable(L, -3);
      lua_pop(L, 1);
    }
  }
  lua_pop(L,1);
} /*}}}*/

void custom_set_node(lua_State* L, Node* node) { /*{{{*/
  int j,i,sz1,sz2;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  if (node->Ndims == 1) {
    for (j=0; j<node->dim; j++) {
      lua_pushnumber(L, j+1);
      lua_pushnumber(L, node->value[j]);
      lua_settable(L, -3);
    }
  } else {
    sz1 = node->size[0]; sz2 = node->size[1];
    for (j=0; j<sz1; j++) {
      lua_pushnumber(L, j+1);
      lua_gettable(L, -2);
      for (i=0; i<sz2; i++) {
        lua_pushnumber(L, i+1);
        lua_pushnumber(L, node->value[sz1*i+j]);
        lua_settable(L, -3);
      }
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1);
}
/*}}}*/

void custom_set_node_block(lua_State* L, Node* node, Block* block) { /*{{{*/
  int i, j, k, ind, sz;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  if (node->Ndims == 1) {
    for (k=0; k<block->dim; k++) {
      if (block->components[k].node != node) continue;
      ind = block->components[k].index;
      lua_pushnumber(L, ind+1);
      lua_pushnumber(L, block->value[k]);
      lua_settable(L, -3);
    }
  } else {
    sz = node->size[0];
    for (k=0; k<block->dim; k++) {
      if (block->components[k].node != node) continue;
      ind = block->components[k].index;
      i = (ind%sz) + 1;
      j = (ind/sz) + 1;
      lua_pushnumber(L, i);
      lua_gettable(L, -2);
      lua_pushnumber(L, j);
      lua_pushnumber(L, block->value[k]);
      lua_settable(L, -3);
      lua_pop(L, 1);
      if (node->kind == SYMMETRIC_KIND) {
        lua_pushnumber(L, j);
        lua_gettable(L, -2);
        lua_pushnumber(L, i);
        lua_pushnumber(L, block->value[k]);
        lua_settable(L, -3);
        lua_pop(L, 1);
      }
    }    
  }
  lua_pop(L, 1); /* node->lua_value */
} /*}}}*/

int custom_read_node(lua_State* L, Node* node) { /*{{{*/
  int j,i,sz1,sz2;
  double val;
  if (!lua_istable(L, -1) && !lua_isuserdata(L, -1)) return EXIT_FAILURE;
  if (node->Ndims == 1) {
    for (j=0; j<node->dim; j++) {
      lua_pushnumber(L, j+1);
      lua_gettable(L, -2);
      if (node->kind == INTEGER_KIND) {
        val = luaL_checkinteger(L, -1);
      } else {
        val = luaL_checknumber(L,-1);
      }
      node->value[j] = val;
      lua_pop(L, 1);
    }
  } else {
    sz1 = node->size[0]; sz2 = node->size[1];
    for (j=0; j<sz1; j++) {
      lua_pushnumber(L, j+1);
      lua_gettable(L, -2);
      if (!lua_istable(L, -1) && !lua_isuserdata(L, -1)) return EXIT_FAILURE;
      for (i=0; i<sz2; i++) {
        lua_pushnumber(L, i+1);
        lua_gettable(L, -2);
        if (node->kind == INTEGER_KIND) {
          val = luaL_checkinteger(L, -1);
        } else {
          val = luaL_checknumber(L,-1);
        }
        node->value[sz1*i+j] = val;
        lua_pop(L, 1);
      }
      lua_pop(L, 1);
    }
  } 
  return EXIT_SUCCESS;
} /*}}}*/

void custom_create(lua_State* L, Node* node) { /*{{{*/
  lua_getglobal(L, MODEL_VNAME);
  lua_getfield(L, -1, node->vname);
  lua_getfield(L, -1, MODEL_INIT_VNAME);
  lua_getfield(L, -1, "copy");
  if (lua_isfunction(L, -1)) {
    /* If there is a .copy function, use it for deep copy */
    lua_insert(L, -2);
    lua_call(L, 1, 1);
    lua_setglobal(L, node->vname);
  } else {
    lua_pop(L, 1); /* .copy */
    /* Else just copy (shallow!) the initial value. */
    if (!QUIET) 
      fprintf(stderr, "Warning: Making a shallow copy of the initial "
              "value of '%s' for sampling.\n", node->vname);
    lua_setglobal(L, node->vname);
  }
  lua_pop(L, 2); /* MODEL_VNAME.node->vname */
} /*}}}*/
