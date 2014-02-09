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

// Test if boundaries need to be adjusted
void checkPartition(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, double partition_time, param *params)
{
    int i;
    fluid_particle *p;
    float h = params->tunable_params.smoothing_radius;

    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    // Get elapsed time since last partition and set new partition time
    double seconds_self =  partition_time;

    float length = params->tunable_params.node_end_x - params->tunable_params.node_start_x;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Send elapsed time to procs to the left and right
    double node[2] = {seconds_self, length};
    double left[2];
    double right[2];
    int tag = 627;
    // Send number of particles to  right and receive from left
    MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_right, tag, left,2,MPI_DOUBLE,proc_to_left,tag,MPI_COMM_COMPUTE,MPI_STATUS_IGNORE);
    // Send number of particles to left and receive from right
    tag = 895;
    MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_left, tag, right,2,MPI_DOUBLE,proc_to_right,tag,MPI_COMM_COMPUTE,MPI_STATUS_IGNORE);

    // Number of seconds for left/right ranks
    double seconds_left = left[0];
    double seconds_right = right[0];

    // Partition length of left/right ranks
    float length_left = left[1];
    float length_right = right[1];

    // Difference in time since last partitioning between self and left/right ranks
    double diff_left = seconds_self - seconds_left;
    double diff_right = seconds_self - seconds_right;

    // Must use the average time as a base otherwise the partitions will not neccessarily move the same ammount
    double average_left = (seconds_self + seconds_left)/2.0;
    double average_right = (seconds_self + seconds_right)/2.0;

    // Allow a 10% difference in the average time
    double max_diff_left = average_left * 0.1;
    double max_diff_right = average_right * 0.1;

    debug_print("max_diff_left %f, max_diff_right %f, diff_left: %f, diff_right %f\n",max_diff_left,max_diff_right,diff_left,diff_right);

    // Adjust left boundary
    // Ensure partition length is atleast 4*h
    if (rank != 0) // Dont move left most boundary
    {
        if( diff_left > max_diff_left && length > 4*h) 
            params->tunable_params.node_start_x += h;
        else if (diff_left < -max_diff_left && length_left > 4*h)
            params->tunable_params.node_start_x -= h;
    }
    // Adjust right boundary
    if (rank != (nprocs-1))
    {
        if( diff_right > max_diff_right&& length > 4*h)
            params->tunable_params.node_end_x -= h;
        else if (diff_right < -max_diff_right && length_right > 4*h)
            params->tunable_params.node_end_x += h;
    }

    debug_print("rank %d node_start %f node_end %f \n", rank, params->tunable_params.node_start_x, params->tunable_params.node_end_x);
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
