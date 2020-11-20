/* -*- mode: C; mode: fold -*- */
/*
 * src/vector_type.c
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

#include "vector_type.h"
#include <stdlib.h>

void vector_set_single_value(lua_State* L, const double value, 
                             const Node* node, const int ind) { /*{{{*/
  int i, j, sz;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  if (node->Ndims == 1) {
    lua_pushnumber(L, value);
    lua_rawseti(L, -2, ind+1);
  } else {
    sz = node->size[0];
    i = (ind%sz)+1; j = (ind/sz)+1;
    lua_rawgeti(L, -1, i);
    lua_pushnumber(L, value);
    lua_rawseti(L, -2, j);
    lua_pop(L, 1);
    if (node->kind == SYMMETRIC_KIND && i != j) {
      lua_rawgeti(L, -1, j);
      lua_pushnumber(L, value);
      lua_rawseti(L, -2, i);
      lua_pop(L, 1);
    }
  }
  lua_pop(L,1);
} /*}}}*/

void vector_set_node(lua_State* L, const Node* node) { /*{{{*/
  int j,i,sz1,sz2;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  if (node->Ndims == 1) {
    for (j=0; j<node->dim; j++) {
      lua_pushnumber(L, node->value[j]);
      lua_rawseti(L, -2, j+1);
    }
  } else {
    sz1 = node->size[0]; sz2 = node->size[1];
    for (j=0; j<sz1; j++) {
      lua_rawgeti(L, -1, j+1);
      for (i=0; i<sz2; i++) {
        lua_pushnumber(L, node->value[sz1*i+j]);
        lua_rawseti(L, -2, i+1);
      }
      lua_pop(L, 1);
    }
  }
  lua_pop(L, 1); /* node->lua_value */
} /*}}}*/

void vector_set_node_block(lua_State* L, const Node* node, 
                           const Block* block) { /*{{{*/
  int i, j, k, ind, sz;
  lua_rawgeti(L, LUA_REGISTRYINDEX, node->lua_value);
  if (node->Ndims == 1) {
    for (k=0; k<block->dim; k++) {
      if (block->components[k].node != node) continue;
      ind = block->components[k].index;
      lua_pushnumber(L, block->value[k]);
      lua_rawseti(L, -2, ind+1);
    }
  } else {
    sz = node->size[0];
    for (k=0; k<block->dim; k++) {
      if (block->components[k].node != node) continue;
      ind = block->components[k].index;
      i = (ind%sz)+1; j = (ind/sz)+1;
      lua_rawgeti(L, -1, i);
      lua_pushnumber(L, block->value[k]);
      lua_rawseti(L, -2, j);
      lua_pop(L, 1);
      if (node->kind == SYMMETRIC_KIND && i != j) {
        lua_rawgeti(L, -1, j);
        lua_pushnumber(L, block->value[k]);
        lua_rawseti(L, -2, i);
        lua_pop(L, 1);
      }
    }    
  }
  lua_pop(L, 1); /* node->lua_value */
} /*}}}*/

int vector_read_node(lua_State* L, Node* node) { /*{{{*/
  int j,i,sz1,sz2;
  double val;
  if (!lua_istable(L,-1)) return EXIT_FAILURE;
  if (node->Ndims == 1) {
    if (lua_objlen(L,-1) != node->dim) return EXIT_FAILURE;
    for (j=0; j<node->dim; j++) {
      lua_rawgeti(L, -1, j+1);
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
    if (lua_objlen(L, -1) != sz1) return EXIT_FAILURE;
    for (j=0; j<sz1; j++) {
      lua_rawgeti(L, -1, j+1);
      if (!lua_istable(L, -1) || lua_objlen(L, -1) != sz2) return -1;
      for (i=0; i<sz2; i++) {
        lua_rawgeti(L, -1, i+1);
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

void vector_create(lua_State* L, Node* node) { /*{{{*/
  int j, i, sz1, sz2;
  if (node->Ndims == 1) {
    lua_createtable(L, node->dim, 0);
    for (j=0; j<node->dim; j++) {
      lua_pushnumber(L, node->value[j]);
      lua_rawseti(L, -2, j+1);
    }
  } else {
    sz1 = node->size[0]; sz2 = node->size[1];
    lua_createtable(L, sz1, 0);
    for (j=0; j<sz1; j++) {
      lua_createtable(L, sz2, 0);
      for (i=0; i<sz2; i++) {
        lua_pushnumber(L, node->value[i*sz1+j]);
        lua_rawseti(L, -2, i+1);
      }
      lua_rawseti(L, -2, j+1);
    }
  }
  lua_setglobal(L, node->vname);
} /*}}}*/
