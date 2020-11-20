/* -*- mode: C; mode: fold -*- */
/*
 * src/matrix_type.h
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

#ifndef _MATRIX_TYPE_H
#define _MATRIX_TYPE_H

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include "grapham_types.h"

/* The functions matrix_* correspond to the functions
 * vector_* in vector_type.h. */

/* Set the value of a single element. */
void matrix_set_single_value(lua_State* L, double value, 
                             Node* node, const int ind);

/* Set the values of the node in given block. */
void matrix_set_node_block(lua_State* L, Node* node, Block* block);

/* Set the value of the whole node. */
void matrix_set_node(lua_State* L, Node* node);

/* Read a value (on top of the stack) to node. */
int matrix_read_node(lua_State* L, Node* node);

/* Create a global value for the node. */
void matrix_create(lua_State* L, Node* node);

#endif
