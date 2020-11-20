/* -*- mode: C; mode: fold -*- */
/*
 * src/number_type.h
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

#ifndef _NUMBER_TYPE_H
#define _NUMBER_TYPE_H

/* Set the Lua value of a variable.
 * 
 * In:
 *   L     -- Lua state.
 *   value -- The value.
 *   node  -- The node. Field used:
 *     vname -- The variable name.
 */
#define number_set_single_value(L,value,node) \
  lua_pushnumber(L, value); \
  lua_setglobal(L, node->vname);

/* Create a Lua value for a variable.
 * 
 * In:
 *   L     -- Lua state.
 *   node  -- The node. Field used:
 *     value -- The value of the variable.
 *     vname -- The variable name.
 */
#define number_create(L, node) \
  lua_pushnumber(L, node->value[0]); \
  lua_setglobal(L, node->vname);

#endif
