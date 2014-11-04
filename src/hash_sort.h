/*
The MIT License (MIT)

Copyright (c) 2014 Adam Simpson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#ifndef fluid_hash_h
#define fluid_hash_h

#include "structs.h"

uint hash_val(float x, float y, float z, neighbor_grid_t *grid, param_t *params);
void hash_particles(fluid_sim_t *fluid_sim);
void sort_hash(fluid_sim_t *fluid_sim);
void find_cell_bounds(fluid_sim_t *fluid_sim);
void fill_particle_neighbors(fluid_sim_t *fluid_sim, uint p_index);
void fill_neighbors(fluid_sim_t *fluid_sim);

// This gets compiled by both c and c++ compilers
#ifdef __cplusplus
extern "C" {
#endif
void find_all_neighbors(fluid_sim_t *fluid_sim);
#ifdef __cplusplus
}
#endif

#endif
