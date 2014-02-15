#include <stdio.h>
#include "geometry.h"
#include "fluid.h"

void constructFluidVolume(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, AABB* fluid, int start_x, 
			  int number_particles_x, edge *edges, float spacing, param *params)
{
    int num_y;
    
    // Number of particles in y,z, number in x is passed in
    num_y = floor((fluid->max_y - fluid->min_y ) / spacing);
    
    // zero out number of edge particles
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;
    
    // Place particles inside bounding volume
    float x,y;
    int nx,ny;
    int i = 0;
    fluid_particle *p;
    for(ny=0; ny<num_y; ny++) {
        y = fluid->min_y + ny*spacing;
        for(nx=0; nx<number_particles_x; nx++) {
            x = fluid->min_x + (start_x + nx)*spacing;
            p = fluid_particles + i;
            p->x = x;
            p->y = y;
            
            // Set pointer array
            fluid_particle_pointers[i] = p;
	    fluid_particle_pointers[i]->id = i;
            i++;
        }
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    printf("rank %d max fluid x: %f\n", rank,fluid->min_x + (start_x + nx-1)*spacing);

    params->number_fluid_particles_local = i;
    params->max_fluid_particle_index = i - 1;
}

// Sets upper bound on number of particles, used for memory allocation
// These numbers are set judiciously for TitanTitan as the number of particles is always small
void setParticleNumbers(AABB *boundary_global, AABB *fluid_global, edge *edges, oob *out_of_bounds, int number_particles_x, float spacing, param *params)
{
    int num_x, num_y, max_y;
    

    // Set fluid local
    num_x = number_particles_x;
    num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    max_y = floor((boundary_global->max_y - boundary_global->min_y ) / spacing);

    // Maximum edge(halo) particles
    edges->max_edge_particles = params->number_fluid_particles_global;

    // The out of bounds particles can become quite large
    // If a flood of particles flows into and then out of a node
    // This will be large
    out_of_bounds->max_oob_particles = params->number_fluid_particles_global;

    // Initial fluid particles
    int num_initial = num_x * num_y;
    printf("initial number of particles %d\n", num_initial);

    // Allow space for all particles if neccessary
    int num_local_max = params->number_fluid_particles_global;

    out_of_bounds->number_vacancies = 0;
    
}

// Set local boundary and fluid particle
void partitionProblem(AABB *boundary_global, AABB *fluid_global, int *x_start, int *length_x, float spacing, param *params)
{
    int i;
    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    // number of fluid particles in x direction
    // +1 added for zeroth particle
    int fluid_particles_x = floor((fluid_global->max_x - fluid_global->min_x ) / spacing) + 1;
    
    // number of particles x direction
    int *particle_length_x = malloc(nprocs*sizeof(int));
    
    // Number of particles in x direction assuming equal spacing
    int equal_spacing = floor(fluid_particles_x/nprocs);
    
    // Initialize each node to have equal width
    for (i=0; i<nprocs; i++)
        particle_length_x[i] = equal_spacing;
    
    // Remaining particles from equal division
    int remaining = fluid_particles_x - (equal_spacing * nprocs);
    
    // Add any remaining particles sequantially to left most nodes
    for (i=0; i<nprocs; i++)
        particle_length_x[i] += (i<remaining?1:0);
    
    // Number of particles to left of current node
    int number_to_left = 0;
    for (i=0; i<rank; i++)
        number_to_left+=particle_length_x[i];
       
    // starting position of nodes x particles
    *x_start = number_to_left;
    // Number of particles in x direction for node
    *length_x = particle_length_x[rank];
        
    // Set node partition values
    params->tunable_params.node_start_x = fluid_global->min_x + ((number_to_left-1) * spacing);
    params->tunable_params.node_end_x   = params->tunable_params.node_start_x + (particle_length_x[rank] * spacing);
    
    if (rank == 0)
        params->tunable_params.node_start_x  = boundary_global->min_x;
    if (rank == nprocs-1)
        params->tunable_params.node_end_x   = boundary_global->max_x;

    // Update requested number of particles with actual value used
    int num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    int total_x = 0;
    for(i=0; i<nprocs; i++)
        total_x += particle_length_x[i];
    params->number_fluid_particles_global = total_x * num_y;

    free(particle_length_x);

}

////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
float min(float a, float b){
    float min = a;
    min = b < min ? b : min;
    return min;
}

float max(float a, float b){
    float max = a;
    max = b > max ? b : max;
    return max;
}

int sgn(float x) {
    int val = 0;
    if (x < 0.0)
        val = -1;
    else if (x > 0.0)
        val = 1;

    return val;
}
