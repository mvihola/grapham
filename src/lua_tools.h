/* -*- mode: C; mode: fold -*- */
/*
 * src/lua_tools.h
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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include <lua.h>

#include "grapham_types.h"

/* Generic error function. */
void error (lua_State *L, const char *fmt, ...);

/* Get the type (including dimension, length, size of each dim.)
 * of the object at index 'ind'. */
Var_type lua_gettype(lua_State* L, int ind, int* dim, int* size, 
                     int* Ndims);

/* Read the length (of first dimension) of the Lua object at 
 * index "ind" */
int lua_getlength(lua_State* L, int ind);

/* Read the vector of length 'Nx' at index 'ind' to array 'x'. */
int lua_getvector(lua_State* L, int ind, double* x, int Nx);

/* Read lower-triangular array from Lua to upper-triangular
 * aray in C. */
int lua_gettril(lua_State* L, int ind, double* x, int Nx);
