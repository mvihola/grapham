/* -*- mode: C; mode: fold -*- */
/*
 * src/lua_math.c
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

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <sys/time.h>

#include "grapham_math.h"
#include "grapham_lib.h"
#include "lua_math.h"

/* Special functions */
static int L_lngamma(lua_State* L) { /*{{{*/
  double x = luaL_checknumber(L, 1);
  lua_pushnumber(L, dlgama_(&x));
  return 1;
} /*}}}*/

/* Special functions */
static int L_systime_usec(lua_State* L) { /*{{{*/
  struct timeval tv;
  gettimeofday(&tv, NULL);
  lua_pushnumber(L, ((double)(tv.tv_sec))*1.0e6 
     + (double)(tv.tv_usec));
  return 1;
} /*}}}*/

/* Just register all functions to Lua globals */
int export_lua_math(lua_State* L) {
  int k;
  
  /* Some constants */
  lua_pushnumber(L, GRAPHAM_INF);
  lua_setglobal(L, "INF");

  lua_pushnumber(L, GRAPHAM_NINF);
  lua_setglobal(L, "NINF");

  /* Special functions */
  lua_pushcfunction(L, L_lngamma);
  lua_setglobal(L, "lngamma");

  lua_pushcfunction(L, L_systime_usec);
  lua_setglobal(L, "systime_usec");

  for(k=0;; k++) {
    if (grapham_dlib[k].name == NULL) {
      break;
    } else {
      lua_pushcfunction(L, grapham_dlib[k].L_fun);
      lua_setglobal(L, grapham_dlib[k].name);
    }
  }
  
  return EXIT_SUCCESS;
}
