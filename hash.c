#include "hash.h"
#include "fluid.h"
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <assert.h>

// Uniform grid hash, this prevents having to check duplicates when inserting
// Fabs needed as neighbor search can go out of bounds
unsigned int hash_val(double x, double y, param *params)
{
    const double spacing = params->smoothing_radius;
    // Calculate grid coordinates
    unsigned int grid_x,grid_y;
    grid_x = floor(x/spacing);
    grid_y = floor(y/spacing);

    // If using glboal boundary size this can be static
    int num_x = params->grid_size_x;
    int num_y = params->grid_size_y;

    unsigned int grid_position = (grid_y * num_x + grid_x);

    return grid_position;
}


// Add halo particles to neighbors array
// Each halo particle looks at the 'behind' fluid particles and adds itself to any within h
void hash_halo(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket *hash, param *params)
{
    int index,i,dx,dy,n, grid_x, grid_y;
    double x,y,r;
    fluid_particle *h_p, *p;
    int n_s = params->number_fluid_particles_local; // Start of halo particles
    int n_f = n_s + params->number_halo_particles;  // End of halo particles
    double spacing = params->smoothing_radius;
    double h = params->smoothing_radius;
    neighbor *ne;

    // Loop over each halo particle
    for(i=n_s; i<n_f; i++)
    {
	// Retrieve hash index for halo particle
        h_p = fluid_particle_pointers[i];

	// Calculate coordinates within bucket grid
	grid_x = floor(h_p->x/spacing);
	grid_y = floor(h_p->y/spacing);

        // Check neighbors of current bucket
        // This only checks 'behind' neighbors as 'forward' neighbors are fluid particles
        for (dx=-1; dx<=0; dx++) {
    	    for (dy=-1; dy<=(dx?1:-1); dy++) {

	        // Calculate index of neighbor cell
	        index = (grid_x + dy)*params->grid_size_x + (grid_y + dx);

                // Go through each fluid particle, p, in neighbor point bucket
                for (n=0;n<hash[index].number_fluid;n++) {
                    p = hash[index].fluid_particles[n];
	
		    // Enforce cutoff
                    r = sqrt((h_p->x-p->x)*(h_p->x-p->x) + (h_p->y-p->y)*(h_p->y-p->y));
                    if(r > h)
                        continue;

                    // Get neighbor bucket for particle p and add halo particle to it
                     ne = &neighbors[p->id];
                     if (ne->number_fluid_neighbors < 300) {
                         ne->fluid_neighbors[ne->number_fluid_neighbors] = h_p;
                         ne->number_fluid_neighbors++;
                     }

                }

            } // End neighbor search y
        } // End neighbor search x

    } // End halo particle loop

} 

// The following function will fill the i'th neighbor bucket with the i'th fluid_particle_pointers particle neighbors
// Only the forward half of the neighbors are added as the forces are symmetrized.
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket * hash, param *params)
{
        int i,j,dx,dy,n,c;
        double x,y, px,py;
        const double spacing = params->smoothing_radius;
	const double h = params->smoothing_radius;
        int n_f = params->number_fluid_particles_local;
        fluid_particle *p, *q, *q_neighbor;
        neighbor *ne;
        double r; 
	unsigned int index, neighbor_index;

        // zero out number of particles in bucket
        for (index=0; index<params->length_hash; index++){
            hash[index].number_fluid = 0;
        }
        
        // First pass - insert fluid particles into hash
        for (i=0; i<n_f; i++) {
            p = fluid_particle_pointers[i];

            neighbors[i].number_fluid_neighbors = 0;
            
            index = hash_val(p->x, p->y, params);

            if (hash[index].number_fluid < 60) {
                hash[index].fluid_particles[hash[index].number_fluid] = p;
                hash[index].number_fluid++;
            }
	    else
		printf("first pass overflow\n");
        }

        // Second pass - fill particle neighbors by processing grid of buckets
        for (j=0; j<params->grid_size_y; j++) {
            for(i=0; i<params->grid_size_x; i++) {

            index = (j * params->grid_size_x + i);
	    if(hash[index].number_fluid == 0)
	        continue;

            // Process current buckets own particle interactions
            // This will only add one neighbor entry per force-pair
            for(c=0; c<hash[index].number_fluid; c++) {
                p = hash[index].fluid_particles[c];
                ne = &neighbors[p->id];
                for(n=c+1; n<hash[index].number_fluid; n++) {
                   q = hash[index].fluid_particles[n];
                   // Append q to p's neighbor list
                    r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y));
                    if(r > h)
                        continue;

                   if(ne->number_fluid_neighbors <300)
                   ne->fluid_neighbors[ne->number_fluid_neighbors++] = q;
                   else
                      printf("self bucket overflow\n");
                }
            }

            // Check neighbors of current bucket
	    // This only checks "forward" neighbors
            for (dx=0; dx<=1; dx++) {
                for (dy=(dx?-1:1); dy<=1; dy++) {
		
	  	    // If the neighbor is outside of the grid we don't process it
		    if ( j+dy < 0 || x+dx < 0 || (i+dx) >= params->grid_size_x || (j+dy) >= params->grid_size_y)
		        continue;

		    neighbor_index = (j+dy)*params->grid_size_x + (i+dx);

			
                    // Add neighbor particles to particles in current bucket
                    for (c=0; c<hash[index].number_fluid; c++) {
		        // Particle in currently being worked on buccket
                        q = hash[index].fluid_particles[c];
                        ne = &neighbors[q->id];
	                for(n=0; n<hash[neighbor_index].number_fluid; n++){
                           // Append neighbor to q's neighbor list
		           q_neighbor = hash[neighbor_index].fluid_particles[n];

                           r = sqrt((q_neighbor->x-q->x)*(q_neighbor->x-q->x) + (q_neighbor->y-q->y)*(q_neighbor->y-q->y));
                           if(r > h)
                               continue;

			    if(ne->number_fluid_neighbors <300)
		            ne->fluid_neighbors[ne->number_fluid_neighbors++] = q_neighbor;
			   else
				printf("neighbor overflow\n");
		        }
                     }
                      
                } // end dy
             }  // end dx
	     
            } // end grid y
        } // end grid x
}// end function
