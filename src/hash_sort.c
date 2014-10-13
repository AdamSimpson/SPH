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

#include "hash_sort.h"
#include <thrust/sort.h>
#include <thrust/execution_policy.h>

// Uniform grid hash, this prevents having to check duplicates when inserting
// Fabs needed if neighbor search gos out of bounds
uint hash_val(float x, float y, param *params)
{
    float spacing = params->grid_spacing;
    float size_x  = params->grid_size_x;

    // Calculate grid coordinates
    uint grid_x,grid_y;
    grid_x = floor(x/spacing);
    grid_y = floor(y/spacing);

    uint grid_position = (grid_y * size_x + grid_x);

    return grid_position;
}

// Hash all particles
void hash_particles(fluid_sim_t *fluid_sim)
{
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
    uint *hash_values = fluid_sim->neighbor_grid->hash_values;
    uint *particle_ids = fluid_sim->neighbor_grid->particle_ids;
    param_t *params = fluid_sim->params;

    int num_particles = params->number_fluid_particles_local + params->number_halo_particles;

    int i;
    for(i=0; i<num_particles; i++) {
        p = fluid_particle_pointers[i];
        hash_values[i] =  hash_val(p->x, p->y, params);
        particle_ids[i] = i;
    }
}

// Sort list of particle id's based upon what their hash value is
void sort_hash(fluid_sim_t *fluid_sim)
{
    param_t *params = fluid_sim->params;
    uint *keys = fluid_sim->neighbor_grid->hash_values;
    uint *values = fluid_sim->neighbor_grid->particle_ids;

    int total_particles = params->number_fluid_particles_local + params->number_halo_particles;

    thrust::sort_by_key(thrust::host, keys, keys+total_particles, values);
}

// Find start and end of hash cells
// Method taken from NVIDIA SDK:
// http://docs.nvidia.com/cuda/samples/5_Simulations/particles/doc/particles.pdf
// Note that end index is one past the "end"
void find_cell_bounds(fluid_sim_t *fluid_sim)
{
    uint *starts = fluid_sim->neighbor_grid->start_indexes;
    uuint *ends = fluid_sim->neighbor_grid->end_indexes;
    uint *hash_values = fluid_sim->neighbor_grid->hash_values;
    uint *particle_ids = fluid_sim->neighbor_grid->particle_ids;
    param_t *params = fluid_sim->params;

    // Reset start indexes
    unsigned int length_hash = params->grid_size_x * params->grid_size_y;
    cudaMemset(starts, 0xffffffff, length_hash*sizeof(uint));

    int num_particles = params->number_fluid_particles_local + params->number_halo_particles;

    int i, i_prev;
    uint hash;
    uint hash_prev;

    // If this particle has a different cell index to the previous
    // particle then it must be the first particle in the cell,
    // so store the index of this particle in the cell.
    // As it isn't the first particle, it must also be the cell end of
    // the previous particle's cell
    for(i=0; i<num_particles; i++) 
    {
        hash = hash_values[i];
        i_prev = (i-1)>=0?(i-1):0;
        hash_prev = hash_values[i_prev]

        if(i==0 || hash!=hash_prev)
        {
            start_indexes[hash] = i;
 
            if(i > 0)
                end_indexes[hash_prev] = i;
        }

        if(i == num_particles - 1)
            end_indexes[hash] = i+1;
    }
}

void fill_particle_neighbors(fluid_sim_t *fluid_sim, fluid_particle_t *p)
{
    // Get neighbor bucket for particle p
    neighbor_t * neighbors = fluid_sim->neighbor_grid->neighbors[p->id];
    uint max_neighbors = fluid_sim->neighbor_grid->max_neighbors;
    param_t *params = fluid_sim->params;

    float smoothing_radius2 = params->tunable_params.smoothing_radius * params->tunable_params.smoothing_radius;

    float spacing = params->grid_spacing;

    // Calculate coordinates within bucket grid
    grid_x = floor(p->x/spacing);
    grid_y = floor(p->y/spacing);

    // Go through neighboring grid buckets
    for(int dy=-1; dy<=1; dy++) {
        for(int dx=-1; dx<=1; dx++) {

            // If the neighbor grid bucket is outside of the grid we don't process it
            if ( grid_y+dy < 0 || grid_x+dx < 0 || (grid_x+dx) >= params->grid_size_x || (grid_y+dy) >= params->grid_size_y)
                continue;

             // Linear hash index for grid bucket
             bucket_index = (grid_y+dy) *params->grid_size_x + grid_x+dx;

             // Start index for hash value of current neighbor grid bucket
             start_index = start_indexes[bucket_index];

             // If neighbor grid bucket is not empty
             if (start_index != 0xffffffff)
             {
                end_index = end_indexes[bucket_index];

                for(int j=start_index; j<end_index; j++)
                {
                    q = fluid_particle_pointers[particle_ids[j]];

                    // Continue if same particle
                    if (p==q)
                        continue;
                    
                    // Calculate distance squared
                    r2 = (p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y);

                    // If inside smoothing radius and enough space in p's neighbor bucket add q
                    if(r2<smoothing_radius2 && neighbors->number_fluid_neighbors < max_neighbors)
                        neighbors->fluid_neighbors[neighbors->number_fluid_neighbors++] = q;
                }
           }
       }
   }

}

void fill_neighbors(fluid_sim_t *fluid_sim)
{

}
