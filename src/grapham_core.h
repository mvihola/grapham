/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_core.h
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
#include "grapham_types.h"

/* Pre-initialisation. L must be OPEN Lua state (that is, you must
 * call  L = lua_open();  before calling this function. */
void grapham_preinit(lua_State *L);

/* After executing the Lua files, this function reads in all the
 * information from the Lua structures, does memory allocations etc. */
void grapham_init(lua_State* L, Model *model, Functional* functional, 
                  Parameters* para);

/* The actual run */
void grapham_run(lua_State* L, Model *model, Functional* functional, 
                  Parameters* para);

/* Some closing routines to be done before quitting. */
void grapham_close(lua_State* L, Model *model, Parameters* para);
