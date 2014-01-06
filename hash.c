#include "hash.h"
#include "fluid.h"
#include <math.h>
#include <stdbool.h>

////////////////////////////////////////////////////////////////////////////
// Spatial Hash
// http://graphics.ethz.ch/Downloads/Publications/Papers/2003/Tes03/Tes03.pdf
////////////////////////////////////////////////////////////////////////////
unsigned int hash_val(double x, double y, double z, double spacing, int hash_size)
{
    unsigned int x_d = (int)floor(x/spacing);
    unsigned int y_d = (int)floor(y/spacing);
    unsigned int z_d = (int)floor(z/spacing);
    const static unsigned int p1 = 73856093;
    const static unsigned int p2 = 19349663;
    const static unsigned int p3 = 83492791;
    
    unsigned int val = (x_d*p1 ^ y_d*p2 ^ z_d*p3) % hash_size;
    
    return val;
}

// Add halo particles to neighbors array
// Halo particles are not added to hash, only neighbors list
// Neighbors may be more than h away...since distance is computed in all smoothing functions
// it is a waste to check as we hash as well
void hash_halo(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket *hash, param *params)
{

    int index,i,dx,dy,dz,dupes,n;
    double x,y,z;
    bool duped;
    fluid_particle *h_p, *q;
    int n_s = params->number_fluid_particles_local;
    int n_f = n_s + params->number_halo_particles;
    double spacing = params->smoothing_radius;
    neighbor *ne;

    for(i=n_s; i<n_f; i++)
    {
        h_p = fluid_particle_pointers[i];
        index = hash_val(h_p->x,h_p->y,h_p->z,spacing,params->length_hash);

        // This is an ugly mess of for loops...
        // This will only find the "lower left" neighbors as forces are symetric
        for (dx=-1; dx<=0; dx++) {
            x = h_p->x + dx*spacing;
            for (dy=-1; dy<=(dx?1:0); dy++) {
                y = h_p->y + dy*spacing;
                for (dz=-1; dz<=((dx|dy)?1:0); dz++) {
                    z = h_p->z + dz*spacing;
                    // Calculate hash index at neighbor point
                    index = hash_val(x,y,z,spacing,params->length_hash);
                      // Go through each fluid particle in neighbor point bucket
                      for (n=0;n<hash[index].number_fluid;n++) {
                          q = hash[index].fluid_particles[n];
                          // Get neighbor ne for particle q
                          ne = &neighbors[q->id];
                          // Make sure not to add duplicate neighbors
                          duped = false;
                          for (dupes=0; dupes < ne->number_fluid_neighbors; dupes++) {
                                if (ne->fluid_neighbors[dupes] == h_p) {
                                    duped = true;
                                    break;
                                }
                            }
                            if (!duped && ne->number_fluid_neighbors < 200) {
                                ne->fluid_neighbors[ne->number_fluid_neighbors] = h_p;
                                ne->number_fluid_neighbors++;
                            }

                      }
                }
            }
        }
    }

} 

// Fill fluid particles into hash
// Neighbors may be more than h away...since distance is computed in all smoothing functions
// it is a waste to check as we hash as well
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket * hash, param *params)
{
    
        int i,dx,dy,dz,n,index,dupes;
        bool duped;
        double x,y,z;
        double spacing = params->smoothing_radius;
	double h = params->smoothing_radius;
        int n_f = params->number_fluid_particles_local;
        fluid_particle *p, *q;
        neighbor *ne;
        double r;        

        // zero out number of particles in bucket
        for (index=0; index<params->length_hash; index++){
            hash[index].number_fluid = 0;
	}
        
        // First pass - insert fluid particles into hash
        for (i=0; i<n_f; i++) {
            p = fluid_particle_pointers[i];
            
            neighbors[i].number_fluid_neighbors = 0;
            
            index = hash_val(p->x,p->y,p->z,spacing,params->length_hash);

            if (hash[index].number_fluid < 200) {
                hash[index].fluid_particles[hash[index].number_fluid] = p;
                hash[index].number_fluid++;
            }
        }

        // Second pass - fill particle neighbors
        for (i=0; i<n_f; i++) {
            p = fluid_particle_pointers[i];
            ne = &neighbors[i];

            // Check in grid around currently particle position for neighbors
	    // This will only find the "right" set of neighbors as forces are symetric
            // This is an ugly mess of for loops...perhaps unroll by hand?
            for (dx=0; dx<=1; dx++) {
                x = p->x + dx*spacing;
                for (dy=(dx?-1:0); dy<=1; dy++) {
                    y = p->y + dy*spacing;
                    for (dz=((dx|dy)?-1:0); dz<=1; dz++) {
                        z = p->z + dz*spacing;
                        // Calculate hash index at neighbor point
                        index = hash_val(x,y,z,spacing,params->length_hash);

                        // Go through each fluid particle in neighbor point bucket
                        for (n=0;n<hash[index].number_fluid;n++) {
                            q = hash[index].fluid_particles[n];
		
			    // This should not be in hash unless pressure/density calculated here as well
                            r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y) + (p->z-q->z)*(p->z-q->z));

			    if(p==q || r > h) // Don't add self to neighbors
                                continue;
 
                            // Make sure not to add duplicate neighbors
		            // If hash collisions occur duplicates may be introduced
			    // This is probably ok to try to remove for quicker simulation
                            duped = false;
                            for (dupes=0; dupes < ne->number_fluid_neighbors; dupes++) {
                               if (ne->fluid_neighbors[dupes] == q) {
                                   duped = true;
                                    break;
                                }
                            }
                            if (!duped) {
                                ne->fluid_neighbors[ne->number_fluid_neighbors] = q;
                                ne->number_fluid_neighbors++;
                            }
                       }
                      
                    }
                }
            }

        }
        
}
