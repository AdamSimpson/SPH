#include <stdio.h>
#include "geometry.h"

void constructFluidVolume(fluid_particle_t *fluid_particles, AABB_t* fluid,
                          int start_x, int number_particles_x, edge_t *edges, param_t *params)
{
    double spacing;
    int num_y;
    int num_z;

    spacing = params->smoothing_radius/2.0;
    // Number of particles in y,z, number in x is passed in
    num_y = floor((fluid->max_y - fluid->min_y ) / spacing);
    num_z = floor((fluid->max_z - fluid->min_z ) / spacing);

    // zero out number of edge particles
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;

    // Place particles inside bounding volume
    double x,y,z;
    int nx,ny,nz;
    int i = 0;
    fluid_particle_t *p;
    for(nz=0; nz<num_z; nz++) {
        z = fluid->min_z + nz*spacing;
        for(ny=0; ny<num_y; ny++) {
            y = fluid->min_y + ny*spacing;
            for(nx=0; nx<number_particles_x; nx++) {
                x = fluid->min_x + (start_x + nx)*spacing;
                p = &fluid_particles[i];
                p->x = x;
                p->y = y;
                p->z = z;
                p->id = i;
                i++;
            }
        }
    }
    params->number_fluid_particles_local = i;
    printf("rank %d max fluid x: %f\n", params->rank,fluid->min_x + (start_x + nx-1)*spacing);
}

// Sets upper bound on number of particles, used for memory allocation
void setParticleNumbers(AABB_t *boundary_global, AABB_t *fluid_global, edge_t *edges,
                        oob_t *out_of_bounds, int number_particles_x, param_t *params)
{
    int num_x;
    int num_y;
    int num_z;

    double spacing = params->smoothing_radius/2.0;

    // Set fluid local
    num_x = number_particles_x;
    num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    num_z = floor((fluid_global->max_z - fluid_global->min_z ) / spacing);

    // Maximum edge particles is a set to 4 particle width y,z slab
    edges->max_edge_particles = 4 * num_y*num_z;
    out_of_bounds->max_oob_particles = 4 * num_y*num_z;

    // Initial fluid particles
    int num_initial = num_x * num_y * num_z;
    printf("initial number of particles %d\n", num_initial);
    int num_extra = num_initial/5;

    // Add initial space, extra space for particle transfers, and left/right out of boudns/halo particles
    params->max_fluid_particles_local = num_initial + num_extra + 2*out_of_bounds->max_oob_particles + 2*edges->max_edge_particles;
    printf("Max fluid particles local: %d\n", params->max_fluid_particles_local);
}

// Set local boundary and fluid particle
void partitionProblem(AABB_t *boundary_global, AABB_t *fluid_global, int *x_start,
                      int *length_x, param_t *params)
{
    int i;
    int nprocs = params->nprocs;
    int rank = params->rank;
    double spacing = params->smoothing_radius/2.0;

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
    params->node_start_x = fluid_global->min_x + ((number_to_left-1) * spacing);
    params->node_end_x   = params->node_start_x + (particle_length_x[rank] * spacing);

    if (rank == 0)
        params->node_start_x  = boundary_global->min_x;
    if (rank == nprocs-1)
        params->node_end_x   = boundary_global->max_x;

    free(particle_length_x);

    printf("rank %d, h %f, x_start %d, num_x %d, start_x %f, end_x: %f\n", rank, spacing, *x_start, *length_x, params->node_start_x, params->node_end_x);

}

// Test if boundaries need to be adjusted
void checkPartition(fluid_particle_t *fluid_particles, oob_t *out_of_bounds, param_t *params)
{

    int i;
    fluid_particle_t *p;
    int num_rank = params->number_fluid_particles_local;
    int rank = params->rank;
    int nprocs = params->nprocs;
    double h = params->smoothing_radius;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of particles and partition length  from right and left
    double length = params->node_end_x - params->node_start_x;
    double node[2]  = {(double)num_rank, length};
    double left[2]  = {0.0, 0.0};
    double right[2] = {0.0, 0.0};
    int tag = 627;
    // Send number of particles to  right and receive from left
    MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_right, tag, left,2,MPI_DOUBLE,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number of particles to left and receive from right
    tag = 895;
    MPI_Sendrecv(node, 2, MPI_DOUBLE, proc_to_left, tag, right,2,MPI_DOUBLE,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Number of particles in left/right ranks
    int num_left = (int)left[0];
    int num_right = (int)right[0];

    // Partition length of left/right ranks
    double length_left = left[1];
    double length_right = right[1];

    int even_particles = params->number_fluid_particles_global/(double)params->nprocs;
    int max_diff = even_particles/10.0f;

    // Difference in particle numbers from an even distribution
    int diff_left  = num_left  - even_particles;
    int diff_right = num_right - even_particles;
    int diff_self  = num_rank  - even_particles;

    // Look at "bins" formed by node start/ends from right to left
    // Only modify node start based upon bin to the left
    // Node end may be modified by node to the rights start
    // Particles per proc if evenly divided

    // current rank has too many particles
    if( diff_self > max_diff && length > 2*h && rank != 0)
        params->node_start_x += h;
    // current rank has too few particles
    else if (diff_self < -max_diff && length_left > 2*h && rank != 0)
        params->node_start_x -= h;

    // Rank to right has too many particles and with move its start to left
    if( diff_right > max_diff && length_right > 2*h && rank != nprocs-1)
        params->node_end_x += h;
    // Rank to right has too few particles and will move its start to right
    else if (diff_right < -max_diff && length > 2*h && rank != nprocs-1)
        params->node_end_x -= h;

    printf("rank %d node_start %f node_end %f \n", rank, params->node_start_x, params->node_end_x);
}

////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
double min(double a, double b){
    double min = a;
    min = b < min ? b : min;
    return min;
}

double max(double a, double b){
    double max = a;
    max = b > max ? b : max;
    return max;
}

int sgn(double x) {
    int val = 0;
    if (x < 0.0)
        val = -1;
    else if (x > 0.0)
        val = 1;

    return val;
}
