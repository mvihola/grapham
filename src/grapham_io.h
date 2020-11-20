/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham_io.h
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
#include <stdio.h>
#include "grapham_types.h"

/* Write the header line with variable names. */
void prepare_outfile(lua_State* L, Model* model, Output* output);

/* Write the header line with adaptation variable names. */
void prepare_outadapt(lua_State* L, Model* model, Parameters* para, Output* output);

/* Read a numeric field from a table (if exists) */
double read_numeric_field(lua_State *L, const char* sname, 
                          const char* fname, double default_value);

/* Read the model structure */
void read_model(lua_State *L, Model* model, Parameters* para);

/* Read the parameters structure */
void read_parameters(lua_State *L, Model *model, Parameters* para);

/* Read the functionals structure */
void read_functional(lua_State *L, Model* model, Parameters* para, 
                     Functional* functional);

/* Write the adaptation parameters (Cholesky factor & scaling) for each block,
 * as well as the block structure, to an output file. */
void write_outcfg(FILE* ofile, Model* model);

/* Graph checking algorithms */
int check_dag(Model* model);
void check_connectivity(Model* model);

/* Just go through the predefined blocks, check sanity
 * and mark each variable 'in_block' flag, if it is
 * already in a block. */
void mark_blocks_used(lua_State *L, Model* model);

/* Add a block for each uninstantiated variable that is
 * not currently attached to any block */
void fill_blocks(lua_State *L, Model* model);

/* Add all the rest free variables to a single block */
void fill_single_block(lua_State *L, Model* model);

/* Add each component not instantiated or already in a sampling block
 * to a separate block */
void fill_single_components(lua_State *L, Model* model);

/* Set up the sampling blocks, according to para.block */
void read_blocks(lua_State *L, Model* model);
