/* -*- mode: C; mode: fold -*- */
/*
 * src/grapham.c
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

/* Includes */ /*{{{*/
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

#include <sys/time.h>

#include "grapham_core.h"

#include "lua_tools.h"
#include "grapham_types.h"

/*}}}*/

/* Just the time difference in seconds */
double tv_diff_sec(struct timeval *tv_start, struct timeval *tv_end) { /*{{{*/
  return (double)(tv_end->tv_sec-tv_start->tv_sec) 
     + (double)(tv_end->tv_usec-tv_start->tv_usec)/1e6;
} /*}}}*/

/* Print some statistics */
void print_statistics(Model* model, Functional* functional, 
                      Parameters* para, double td) { /*{{{*/
  long niter = para->niter;
  long nburn = para->nburn;
  long nthin = para->nthin;
  long ntot = niter+nburn;
  long Na1;
  int i,j,k,l, sz, sdim = 0;
  Block* block;
  Node* node;
  boolean full_var;
  
  if (!QUIET) {
    for (l=0; l<model->Nblocks; l++) {
      sdim += model->blocks[l]->dim;
    }
    printf("Total dimension: %i (sampling %i with ", model->dim, sdim);
    switch(para->alg) {
    case ALG_AM: 
      printf("AM"); break;
    case ALG_AMS:
      printf("AMS"); break;
    case ALG_RBAM: 
      printf("RBAM"); break;
    case ALG_RBAMS:
      printf("RBAMS"); break;
    case ALG_ASCM:
      printf("ASCM"); break;
    case ALG_ASHM:
      printf("ASHM"); break;
    case ALG_METROPOLIS:
      printf("Metropolis"); break;
    }
    if (!LUA_CALLS) {
      printf("); No Lua calls\n");
    } else {
      printf("); With Lua calls\n");
    }
    if (nthin>1) {
      printf("Results after %lu samples (%lu total simulated, %lu burn-in):\n",
             functional->called, niter, nburn);
    } else {
      printf("Results after %lu samples (%lu burn-in):\n",niter,nburn);
    }
  }
  
  if (functional->evaluate) {
  if (!QUIET) {
    printf("Functional average = [");
    for (i=0; i<functional->mean_N; i++) printf(" %f",functional->mean[i]);
    printf(" ]\n"); 
  } else {
    if (functional->mean_N>=1) {
      printf(OUTFILE_NUMFORMAT,functional->mean[0]);
      for (i=1; i<functional->mean_N; i++) 
        printf(OUTFILE_FIELD_SEP OUTFILE_NUMFORMAT,functional->mean[i]);
      printf("\n");
    }
  }
  }
  
  if (!QUIET) {
    if (model->Nblocks == 1) {
      block = model->blocks[0];
      if (block->accepted+block->rejected == ntot) {
      printf("Acceptance rate: %.2f%%\n", 
             (double)(block->accepted)/((double)(ntot))*100 );
      } else {
        /* We have delayed rejections. */
        Na1 = 2*ntot-(block->accepted+block->rejected);
        printf("Acceptance rate: %.2f%% (%.2f%%/%.2f%%)\n", 
               (double)(block->accepted)/((double)(ntot))*100,
               (double)(Na1)/((double)(ntot))*100,
               (double)(block->accepted-Na1)/((double)(ntot))*100
               );
      }
    } else {
      printf("Acceptance rates:");
      
      for (l=0; l<model->Nblocks; l++) {
        block = model->blocks[l];
        printf("  ( ");
        for (j=0; j<block->Nvars; j++) {
          node = block->variables[j];
          full_var = TRUE;
          /* Check whether all components of the variable 
           * are in this block */
          for (k=0; k<node->dim; k++) {
            for (i=0; i<block->dim; i++) {
              if (block->components[i].node == node
                  && block->components[i].index == k)
              goto next_elem;
            }
            full_var = FALSE;
            next_elem:
            ;
          }
          if (full_var) {
            /* If all components are in this block: */
            printf("%s ", node->vname);
          } else {
            /* The variable is splitted among blocks. */
            for (i=0; i<block->dim; i++) {
              if (block->components[i].node != node) continue;
              if (node->Ndims == 1) {
                printf("%s[%i] ", node->vname, block->components[i].index+1);
              } else {
                sz = node->size[0];
                printf("%s[%i][%i] ", node->vname, 
                       block->components[i].index%sz+1,
                       block->components[i].index/sz+1);
              }
            }
          }
        }
        printf("): ");
        if (block->accepted+block->rejected == ntot) {
          printf("%.2f%%", 
                 (double)(block->accepted)/((double)(ntot))*100 );
        } else {
          /* We have delayed rejections. */
          Na1 = 2*ntot-(block->accepted+block->rejected);
          printf("%.2f%% (%.2f%%/%.2f%%)", 
                 (double)(block->accepted)/((double)(ntot))*100,
                 (double)(Na1)/((double)(ntot))*100,
                 (double)(block->accepted-Na1)/((double)(ntot))*100
                 );
        }
      }
      printf("\n");
    }
  
    fprintf(stdout, "Sampling CPU time: %.2fs  (%.2fus/sample; %.2fus/sample/dim)\n",
            td, 1e6*td/( (double)(ntot) ), 
            1e6*td/( (double)(ntot) )/((double)sdim));
  }
}

/*}}}*/

int main(int argc, char *argv[]) { /*{{{*/
  int f, lua_inline = -1;
  struct timeval tv_start, tv_end;
  boolean no_lua_input = TRUE;
  Model model;
  Functional functional;
  Parameters para;
  lua_State *L = NULL;

  VERBOSE = FALSE; MOREVERBOSE = FALSE; QUIET = FALSE;
  /* For now... */
  LUA_CALLS = FALSE;
  
  
  /* Handle simple parameters and check for consistency */
  for (f=1; f<argc; f++) {
    if (strcmp(argv[f],"-v") == 0) {
      VERBOSE = TRUE;
    } else if (strcmp(argv[f],"-vv") == 0) {
      VERBOSE = TRUE; 
      MOREVERBOSE = TRUE;
    } else if (strcmp(argv[f],"-q") == 0) {
      VERBOSE = FALSE; 
      MOREVERBOSE = FALSE;
      QUIET = TRUE;
    } else if (strcmp(argv[f],"-e") == 0) {
      if (f>=argc-1) continue;
      no_lua_input = FALSE;
    } else {
      no_lua_input = FALSE;
    }
  }
  if (!QUIET) {
    fprintf(stderr, VERSION_TEXT);
  }
  if (no_lua_input) {
    fprintf(stderr, USAGE_TEXT, argv[0]);
    return EXIT_FAILURE;
  }
  
  /* Pre-initialisation: open standard libraries, register functions etc. */
  L = lua_open();
  grapham_preinit(L);
  
  /* Run the Lua files (and inlines) in the given order. */
  for (f=1; f<argc; f++) {
    if (strcmp(argv[f],"-v") == 0 
        || strcmp(argv[f],"-vv") == 0
        || strcmp(argv[f],"-q") == 0) {
      continue;
    } else if (strcmp(argv[f],"-e") == 0) {
      if (f>=argc-1) {
        fprintf(stderr, USAGE_TEXT, argv[0]);
        return EXIT_FAILURE;
      } else {
        lua_inline = ++f;
        /* Load inline */
        if (VERBOSE) fprintf(stderr, "Executing Lua inline\n");
        if (luaL_loadbuffer(L, argv[lua_inline], 
                            strlen(argv[lua_inline]), "=<inline>") 
            || lua_pcall(L, 0, 0, 0))
        error(L, "Error: %s", lua_tostring(L, -1));
      }
    } else {
      /* Load the model file */
      if (VERBOSE) fprintf(stderr, "Loading the model file...\n");
      if (luaL_loadfile(L, argv[f]) || lua_pcall(L, 0, 0, 0))
        error(L, "Error: %s", lua_tostring(L, -1));
    }
  }

  /* Initialise */
  grapham_init(L, &model, &functional, &para);
  
  /* The actual sampling & estimation */
  gettimeofday(&tv_start, NULL);
  grapham_run(L, &model, &functional, &para);
  gettimeofday(&tv_end, NULL);
  
  /* Do some finalising routines before quitting. */
  grapham_close(L, &model, &para);
  
  /* Some output */
  print_statistics(&model, &functional, &para, 
                   tv_diff_sec(&tv_start, &tv_end));
  
  return EXIT_SUCCESS;
}

/*}}}*/
