/* -*- mode: C; mode: fold -*- */
/*
 * src/vector_type.h
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

#ifndef _VECTOR_TYPE_H
#define _VECTOR_TYPE_H

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include "grapham_types.h"

/* Set the value of a single element.
 * 
 * In:
 *   L     -- Lua state.
 *   value -- The value of the state.
 *   node  -- The variable. Fields used:
 *     lua_value -- The index of the Lua value.
 *     Ndims     -- Number of dimensions (1 or 2).
 *     size      -- The size of each dimension.
 *     kind      -- The kind of the variable.
 *   ind   -- The index of the element to be set.
 */
void vector_set_single_value(lua_State* L, const double value, 
                             const Node* node, const int ind);

/* Set the value of the whole node.
 * 
 * In:
 *   L     -- Lua state.
 *   value -- The value of the state.
 *   node  -- The variable. Fields used:
 *     lua_value -- The index of the Lua value.
 *     dim       -- The length of the node.
 *     Ndims     -- Number of dimensions (1 or 2).
 *     size      -- The size of each dimension.
 *   ind   -- The index of the element to be set.
 */
void vector_set_node(lua_State* L, const Node* node);

/* Set the values of the node in given block.
 * 
 * In:
 *   L     -- Lua state.
 *   node  -- The variable. Fields used:
 *     lua_value -- The index of the Lua value.
 *     Ndims     -- Number of dimensions (1 or 2).
 *     size      -- The size of each dimension.
 *     kind      -- The kind of the variable.
 *   block -- The block whose values are updated to the node.
 *     dim        -- The dimension of the block.
 *     value      -- The value of the block.
 *     components.node, components.index -- The node and the
 *                   index of the corresponding component.
 */
void vector_set_node_block(lua_State* L, const Node* node, const Block* block);

/* Read a value (on top of the stack) of a node.
 * 
 * In:
 *   L     -- The Lua state. It is assumed that the value is
 *            on top of the Lua stack!
 *   node  -- The variable. Fields used:
 *     dim       -- The length of the node.
 *     Ndims     -- Number of dimensions (1 or 2).
 *     size      -- The size of each dimension.
 *     kind      -- The kind of the variable.
 * 
 * Out:
 *   node.value  -- The read value of the node.
 */
int vector_read_node(lua_State* L, Node* node);

/* Create a global value corresponding to a node.
 * 
 * In:
 *   L     -- The Lua state.
 *   node  -- The variable. Fields used:
 *     dim       -- The length of the node.
 *     Ndims     -- Number of dimensions (1 or 2).
 *     value     -- The node value.
 *     vname     -- The Lua variable name.
 */
void vector_create(lua_State* L, Node* node);

#endif
