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

#include "hash.h"
#include "fluid.h"
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <assert.h>

// Uniform grid hash, this prevents having to check duplicates when inserting
// Fabs needed as neighbor search can go out of bounds
unsigned int hash_val(float x, float y, neighbor_grid_t *grid, param *params)
{
    const float spacing = grid->spacing;

    // Calculate grid coordinates
    unsigned int grid_x,grid_y;
    grid_x = floor(x/spacing);
    grid_y = floor(y/spacing);

    unsigned int grid_position = (grid_y * grid->size_x + grid_x);

    return grid_position;
}

// Add halo particles to neighbors array
// We also calculate the density as it's convenient
void hash_halo(fluid_particle **fluid_particle_pointers,  neighbor_grid_t *grid, param *params)
{
    int index,i,dx,dy,n, grid_x, grid_y;
    float x,y,r2, r;
    fluid_particle *h_p, *p;

    int n_start = params->number_fluid_particles_local; // Start of halo particles
    int n_finish = n_start + params->number_halo_particles;  // End of halo particles
    float h = params->tunable_params.smoothing_radius;
    float h_recip = 1.0f/h;

    unsigned int max_neighbors = grid->max_neighbors;
    float spacing = grid->spacing;
    neighbor *neighbors = grid->neighbors;
    bucket_t *grid_buckets = grid->grid_buckets;

    float ratio;
    float h2 = h*h;
    neighbor *ne;

    // Loop over each halo particle
    for(i=n_start; i<n_finish; i++)
    {
	// Retrieve hash index for halo particle
        h_p = fluid_particle_pointers[i];

	// Calculate coordinates within bucket grid
	grid_x = floor(h_p->x_star/spacing);
	grid_y = floor(h_p->y_star/spacing);

        // Check neighbors of current bucket
        // This only checks 'behind' neighbors as 'forward' neighbors are fluid particles
        for (dx=-1; dx<=1; dx++) {
    	    for (dy=-1; dy<=1; dy++) {

                // If the neighbor is outside of the grid we don't process it
                if ( grid_y+dy < 0 || grid_x+dx < 0 || (grid_x+dx) >= grid->size_x || (grid_y+dy) >= grid->size_y)
                    continue;

	        // Calculate index of neighbor cell
	        index = (grid_y + dy)*grid->size_x + (grid_x + dx);

                // Go through each fluid particle, p, in neighbor point bucket
                for (n=0; n<grid_buckets[index].number_fluid; n++) {
                    p = grid_buckets[index].fluid_particles[n];
	
		    // Enforce cutoff
                    r2 = (h_p->x_star-p->x_star)*(h_p->x_star-p->x_star) + (h_p->y_star-p->y_star)*(h_p->y_star-p->y_star);
                    if(r2 > h2)
                        continue;
	
                     // Get neighbor bucket for particle p and add halo particle to it
                     ne = &neighbors[p->id];
                     if (ne->number_fluid_neighbors < max_neighbors) {
                         ne->fluid_neighbors[ne->number_fluid_neighbors++] = h_p;
                     }
		     else
			debug_print("halo overflowing\n");

                }

            } // End neighbor search y
        } // End neighbor search x

    } // End halo particle loop

} 

// The following function will fill the i'th neighbor bucket with the i'th fluid_particle_pointers particle neighbors
// Only the forward half of the neighbors are added as the forces are symmetrized.
// We also calculate the density as it's convenient
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor_grid_t *grid, param *params)
{
        int i,j,dx,dy,n,c;
        float x,y, px,py;
        float h = params->tunable_params.smoothing_radius;
        float h_recip = 1.0f/h;
        float h2 = h*h;
        int n_f = params->number_fluid_particles_local;

        unsigned int max_neighbors = grid->max_neighbors;
        unsigned int max_bucket_size = grid->max_bucket_size;
        neighbor *neighbors = grid->neighbors;
        bucket_t *grid_buckets = grid->grid_buckets; 
        unsigned int length_hash = grid->size_x * grid->size_y;

        fluid_particle *p, *q, *q_neighbor;
        neighbor *ne;
        float r,r2, ratio; 
        unsigned int index, neighbor_index;

        // zero out number of particles in bucket
        for (index=0; index<length_hash; index++){
            grid_buckets[index].number_fluid = 0;
        }
        
        // First pass - insert fluid particles into hash
        for (i=0; i<n_f; i++) {
            p = fluid_particle_pointers[i];

            neighbors[i].number_fluid_neighbors = 0;
            
            index = hash_val(p->x_star, p->y_star, grid, params);

            if (grid_buckets[index].number_fluid < max_bucket_size) {
                grid_buckets[index].fluid_particles[grid_buckets[index].number_fluid] = p;
                grid_buckets[index].number_fluid++;
            }
	    else
		debug_print("first pass overflow\n");
        }

        // Second pass - fill particle neighbors by processing grid of buckets
        for (j=0; j<grid->size_y; j++) {
            for(i=0; i<grid->size_x; i++) {

            index = (j * grid->size_x + i);
	    if(grid_buckets[index].number_fluid == 0)
	        continue;

            // Process current buckets own particle interactions
            // This will only add one neighbor entry per force-pair
            for(c=0; c<grid_buckets[index].number_fluid; c++) {
                p = grid_buckets[index].fluid_particles[c];
                ne = &neighbors[p->id];
                for(n=0; n<grid_buckets[index].number_fluid; n++) {
                   if(c==n)
                       continue;
                   q = grid_buckets[index].fluid_particles[n];
                   // Append q to p's neighbor list
                    r2 = (p->x_star-q->x_star)*(p->x_star-q->x_star) + (p->y_star-q->y_star)*(p->y_star-q->y_star);
                    if(r2 > h2)
                        continue;

                   if(ne->number_fluid_neighbors < max_neighbors) {
                       ne->fluid_neighbors[ne->number_fluid_neighbors++] = q;
                   }
                   else
                      debug_print("self bucket overflow\n");
                }
            }

            // Check neighbors of current bucket
            // This only checks "forward" neighbors
            for (dx=-1; dx<=1; dx++) {
                for (dy=-1; dy<=1; dy++) {
		    if(dx==0 && dy==0)
                        continue;
	  	    // If the neighbor is outside of the grid we don't process it
		    if ( j+dy < 0 || i+dx < 0 || (i+dx) >= grid->size_x || (j+dy) >= grid->size_y)
		        continue;

		    neighbor_index = (j+dy)*grid->size_x + (i+dx);
			
                    // Add neighbor particles to particles in current bucket
                    for (c=0; c<grid_buckets[index].number_fluid; c++) {
		        // Particle in currently being worked on buccket
                        q = grid_buckets[index].fluid_particles[c];
                        ne = &neighbors[q->id];
	                for(n=0; n<grid_buckets[neighbor_index].number_fluid; n++){
                            // Append neighbor to q's neighbor list
		            q_neighbor = grid_buckets[neighbor_index].fluid_particles[n];
                             r2 = (q_neighbor->x-q->x)*(q_neighbor->x-q->x) + (q_neighbor->y-q->y)*(q_neighbor->y-q->y);
                            if(r2 > h2)
                                continue;
                            if(ne->number_fluid_neighbors < max_neighbors) {
		                ne->fluid_neighbors[ne->number_fluid_neighbors++] = q_neighbor;
                             }
                             else
                                 debug_print("neighbor overflow\n");
		        }
                    }
                      
                } // end dy
             }  // end dx
	     
            } // end grid y
        } // end grid x

}// end function
