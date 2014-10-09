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

#include <stdbool.h>

typedef struct BUCKET_T bucket_t;
typedef struct NEIGHBOR_GRID_T neighbor_grid_t;

#include "fluid.h"

struct BUCKET_T {
    fluid_particle **fluid_particles;
    unsigned int number_fluid;
}; // neighbor 'bucket' for hash value

struct NEIGHBOR_GRID_T {
    float spacing;  // Spacing between buckets
    unsigned int size_x; // Number of buckets in x
    unsigned int size_y; // Number of buckets in y
    neighbor *neighbors; // Particle neighbor buckets
    bucket_t *grid_buckets; // Grid to place hashed particles into
    unsigned int max_neighbors; // Maximum neighbors allowed for each particle
    unsigned int max_bucket_size; // Maximum particles in hash bucket
};

unsigned int hash_val(float x, float y, neighbor_grid_t *grid, param *params);
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor_grid_t *grid, param *params);
void hash_halo(fluid_particle **fluid_particle_pointers,  neighbor_grid_t *grid, param *params);

#endif

