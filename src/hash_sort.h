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


typedef struct NEIGHBOR_T neighbor_t;
typedef struct NEIGHBOR_GRID_T neighbor_grid_t;

#include "fluid.h"

// Bucket to hold each particles nearest neighbors
struct NEIGHBOR_T{
    fluid_particle_t **fluid_neighbors;
    int number_fluid_neighbors;
};

struct NEIGHBOR_GRID_T {
    float spacing;  // Spacing between buckets
    uint size_x; // Number of buckets in x
    uint size_y; // Number of buckets in y
    uint *start_indexes; // Start index for hash values
    uint *end_indexes;   // End index for hash values
    uint *hash_values; // Array of hash values
    uint *particle_ids; // Array of particle id's
    uint max_neighbors; // Maximum neighbors allowed for each particle
    neighbor_t *neighbors; // Particle neighbor buckets
};

uint hash_val(float x, float y, neighbor_grid_t *grid, param_t *params);
void hash_particles(fluid_sim_t *fluid_sim);
void sort_hash(fluid_sim_t *fluid_sim);
void find_cell_bounds(fluid_sim_t *fluid_sim);
void fill_particle_neighbors(fluid_sim_t *fluid_sim, fluid_particle_t *p);
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
