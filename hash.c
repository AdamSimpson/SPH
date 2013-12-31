#include "hash.h"
#include "fluid.h"
#include <math.h>
#include <stdint.h>
#include <time.h>

// Uniform grid hash, this prevents having to check duplicates when inserting
uint64_t hash_val(double x, double y, double z, param *params)
{
    double spacing = params->smoothing_radius;
    // Calculate grid coordinates
    uint64_t grid_x,grid_y,grid_z;
    grid_x = floor(fabs(x)/spacing);
    grid_y = floor(fabs(y)/spacing);
    grid_z = floor(fabs(z)/spacing);

    int num_x = params->grid_size_x;
    int num_y = params->grid_size_y;
    int num_z = params->grid_size_z;

    uint64_t grid_position = (num_x * num_y * grid_z) + (grid_y * num_x + grid_x);
    return grid_position;
}

// Asscending hash value comparison for qsort
int compare_hash(const void *a, const void *b)
{
    return ((uint2*)a)->x - ((uint2*)b)->x;
}

void fill_hash_positions(uint2 *hash, uint2 *hash_positions, int hash_length)
{
    uint64_t i;
    uint64_t val,prev_val;

    // create start/end array
    for (i=0; i<hash_length; i++) {
        val = hash[i].x;

        // If the value is not the first and it's different from the previous val it is the start
        if (i == 0 || val != hash[i-1].x) {
            hash_positions[val].x = i;

            // If this is the start of a new hash value the previous one has ended
            if(i > 0) {
                prev_val = hash[i-1].x;
                hash_positions[prev_val].y = i-1;
            }
        }
        // Set end value of last hash value
        if (i == hash_length - 1)
            hash_positions[val].y = i;
    }
}

// Hash and fill fluid particles and place fluid particles into neighbors array
void hash_fluid(fluid_particle* fluid_particles, neighbor *neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, param *params)
{
    uint64_t index, n, hash_start, hash_end, val;
    double x,y,z;
    int i,dx,dy,dz;
    int n_f = params->number_fluid_particles;
    fluid_particle *p, *q;
    boundary_particle *k;
    double spacing = params->smoothing_radius;
    neighbor *ne;

    uint64_t grid_size = (uint64_t)params->grid_size_x * params->grid_size_y * params->grid_size_z;

    // First pass - insert fluid particles into hash
    for (i=0; i<n_f; i++) {
        p = &fluid_particles[i];
        val = hash_val(p->x,p->y,p->z,params);

        // Fill hash: (hash_value, particle ID)
        fluid_hash[i].x = val;
        fluid_hash[i].y = i;

        // Zero out neighbor values
        neighbors[i].number_fluid_neighbors = 0;
    }

    // sort hash by hash value in asscending order
    qsort(fluid_hash, n_f, sizeof(uint2), compare_hash);

    // Invalidate old fluid_hash starting positions
    // Using max is probably not great but should work
    for (i=0; i<grid_size; i++) 
        fluid_hash_positions[i].x = UINT64_MAX;

    // fill fluid hash start/end positions
    fill_hash_positions(fluid_hash, fluid_hash_positions, n_f);

    // Second pass - fill particle neighbors
    for (i=0; i<n_f; i++) {
        p = &fluid_particles[i];
        ne = &neighbors[i];
        // Check in grid around currently particle position for neighbors
        for (dx=-2; dx<=2; dx++) {
            x = p->x + dx*spacing;
            for (dy=-2; dy<=2; dy++) {
                y = p->y + dy*spacing;
                for (dz=-2; dz<=2; dz++) {
                    z = p->z + dz*spacing;
                    // Calculate hash index at neighbor point
                    index = hash_val(x,y,z,params);

                    hash_start = fluid_hash_positions[index].x;
		    hash_end = fluid_hash_positions[index].y;

		    if(hash_start != UINT64_MAX) {
                        // Go through each fluid particle in neighbor point bucket
                        for (n=hash_start;n<=hash_end;n++) {
                            q = &fluid_particles[fluid_hash[n].y];
                            double distance = sqrt((p->x-q->x)*(p->x-q->x)+(p->y-q->y)*(p->y-q->y)+(p->z-q->z)*(p->z-q->z));
                            // Make sure the distance is less than 2h
                            if (distance < 2.0*spacing && ne->number_fluid_neighbors < 30) {
                                ne->fluid_neighbors[ne->number_fluid_neighbors] = q;
                                ne->number_fluid_neighbors++;
                            }
                        }
		    }
                        
                } // dz loop
            } // dy loop
        } // dx loop
    } // particle loop
}
