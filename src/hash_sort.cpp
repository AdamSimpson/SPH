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

// Calculate and fill all neighbor particles
extern "C" void find_all_neighbors(fluid_sim_t *fluid_sim)
{
    hash_particles(fluid_sim);
    sort_hash(fluid_sim);
    find_cell_bounds(fluid_sim);
    fill_neighbors(fluid_sim);
}

// Uniform grid hash, this prevents having to check duplicates when inserting
// Fabs needed if neighbor search gos out of bounds
uint hash_val(float x, float y, float z, fluid_sim_t *fluid_sim)
{
    float spacing = fluid_sim->neighbor_grid->spacing;
    float size_x  = fluid_sim->neighbor_grid->size_x;
    float size_y = fluid_sim->neighbor_grid->size_y;

    // Calculate grid coordinates
    uint grid_x,grid_y,grid_z;
    grid_x = floor(x/spacing);
    grid_y = floor(y/spacing);
    grid_z = floor(z/spacing);

    uint grid_position = (grid_z *(size_x + size_y) + grid_y * size_x + grid_x);

    return grid_position;
}

// Hash all particles
void hash_particles(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    uint *hash_values = fluid_sim->neighbor_grid->hash_values;
    uint *particle_ids = fluid_sim->neighbor_grid->particle_ids;
    param_t *params = fluid_sim->params;

    int num_particles = params->number_fluid_particles_local + params->number_halo_particles;

    uint p_index;

    int i;
    for(i=0; i<num_particles; i++) {
        p_index = fluid_particle_indices[i];
        hash_values[i] =  hash_val(fluid_particles->x[p_index], fluid_particles->y[p_index], fluid_particles->z[p_index], fluid_sim);
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
    uint *start_indices = fluid_sim->neighbor_grid->start_indices;
    uint *end_indices = fluid_sim->neighbor_grid->end_indices;
    uint *hash_values = fluid_sim->neighbor_grid->hash_values;
    uint *particle_ids = fluid_sim->neighbor_grid->particle_ids;
    param_t *params = fluid_sim->params;

    // Reset start indicies
    unsigned int length_hash = fluid_sim->neighbor_grid->size_x 
                             * fluid_sim->neighbor_grid->size_y
                             * fluid_sim->neighbor_grid->size_z;
    memset(start_indices, ((uint)-1), length_hash*sizeof(uint));

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
        hash_prev = hash_values[i_prev];

        if(i==0 || hash!=hash_prev)
        {
            start_indices[hash] = i;
 
            if(i > 0)
                end_indices[hash_prev] = i;
        }

        if(i == num_particles - 1)
            end_indices[hash] = i+1;
    }
}

// Neighbors are accessed multiple times per step so we keep them in buckets
void fill_particle_neighbors(fluid_sim_t *fluid_sim, uint p_index)
{

    uint max_neighbors = fluid_sim->neighbor_grid->max_neighbors;
    param_t *params = fluid_sim->params;
    float smoothing_radius2 = params->tunable_params.smoothing_radius * params->tunable_params.smoothing_radius;

    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;

    uint *particle_ids = fluid_sim->neighbor_grid->particle_ids;

    // Get neighbor bucket for particle p
    neighbor_t * neighbors = &fluid_sim->neighbor_grid->neighbors[fluid_particles->id[p_index]];

    float spacing = fluid_sim->neighbor_grid->spacing;
    uint start_index, end_index;
    int dx,dy,dz, grid_x, grid_y, grid_z, bucket_index;
    uint q_index;

    uint *start_indices = fluid_sim->neighbor_grid->start_indices;
    uint *end_indices = fluid_sim->neighbor_grid->end_indices;

    neighbors->number_fluid_neighbors = 0;

    // Calculate coordinates within bucket grid
    grid_x = floor(fluid_particles->x[p_index]/spacing);
    grid_y = floor(fluid_particles->y[p_index]/spacing);
    grid_z = floor(fluid_particles->z[p_index]/spacing);

    float r2;

    // Go through neighboring grid buckets
    for(dz=-1; dz<=1; dz++) {
        for(dy=-1; dy<=1; dy++) {
            for(dx=-1; dx<=1; dx++) {
                // If the neighbor grid bucket is outside of the grid we don't process it
                if ( grid_y+dy < 0 || grid_x+dx < 0 || grid_z+dz < 0||
                   (grid_x+dx) >= fluid_sim->neighbor_grid->size_x 
                   || (grid_y+dy) >= fluid_sim->neighbor_grid->size_y
                   || (grid_z+dz) >= fluid_sim->neighbor_grid->size_z)
                       continue;

                 // Linear hash index for grid bucket
                 bucket_index = (grid_z+dz)*(fluid_sim->neighbor_grid->size_x + fluid_sim->neighbor_grid->size_y) 
                              + (grid_y+dy) *fluid_sim->neighbor_grid->size_x + grid_x+dx;

                 // Start index for hash value of current neighbor grid bucket
                 start_index = start_indices[bucket_index];

                 // If neighbor grid bucket is not empty
                 if (start_index != (uint)-1)
                 {
                    end_index = end_indices[bucket_index];

                    for(int j=start_index; j<end_index; j++)
                    {
                        q_index = fluid_particle_indices[particle_ids[j]];

                        // Continue if same particle
                        if (p_index==q_index)
                            continue;
                    
                        // Calculate distance squared
                        r2 = (fluid_particles->x[p_index]-fluid_particles->x[q_index])
                            *(fluid_particles->x[p_index]-fluid_particles->x[q_index])
                            +(fluid_particles->y[p_index]-fluid_particles->y[q_index])
                            *(fluid_particles->y[p_index]-fluid_particles->y[q_index])
                            +(fluid_particles->z[p_index]-fluid_particles->z[q_index])
                            *(fluid_particles->z[p_index]-fluid_particles->z[q_index]);

                        // If inside smoothing radius and enough space in p's neighbor bucket add q
                        if(r2<smoothing_radius2 && neighbors->number_fluid_neighbors < max_neighbors)
                            neighbors->fluid_neighbors[neighbors->number_fluid_neighbors++] = q_index;
                    }
               }
            } //dx
       } //dy
   } //dz

}

void fill_neighbors(fluid_sim_t *fluid_sim)
{
    param_t *params = fluid_sim->params;
    int total_particles = params->number_fluid_particles_local;
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;

    int i;
    uint p_index;

    // Fill neighbor bucket for all resident particles
    for(i=0; i<total_particles; i++) {
        p_index = fluid_particle_indices[i];
        fill_particle_neighbors(fluid_sim, p_index);
    }       
}
