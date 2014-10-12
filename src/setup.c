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

#include <stdio.h>
#include "setup.h"
#include "fluid.h"

// Allocate main simulation memory
void alloc_sim(fluid_sim_t *fluid_sim)
{
    unsigned int i;

    // Neighbor grid setup
    fluid_sim->neighbor_grid = (neighbor_grid_t*)malloc(sizeof(neighbor_grid_t));
    fluid_sim->neighbor_grid->max_bucket_size = 100;
    fluid_sim->neighbor_grid->max_neighbors = fluid_sim->neighbor_grid->max_bucket_size*4;
    fluid_sim->neighbor_grid->spacing = fluid_sim->params->tunable_params.smoothing_radius;

    size_t total_bytes = 0;
    size_t bytes;
    // Allocate fluid particles array
    bytes = fluid_sim->params->max_fluid_particles_local * sizeof(fluid_particle_t);
    total_bytes+=bytes;
    fluid_sim->fluid_particles = malloc(bytes);
    if(fluid_sim->fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");

    // Allocate (x,y) coordinate array, transfer pixel coords
    bytes = 2 * fluid_sim->params->max_fluid_particles_local * sizeof(short);
    total_bytes+=bytes;
    fluid_sim->fluid_particle_coords = malloc(bytes);
    if(fluid_sim->fluid_particle_coords == NULL)
        printf("Could not allocate fluid_particle coords\n");

    // Allocate pointer array used to traverse non vacant particles
    bytes = fluid_sim->params->max_fluid_particles_local * sizeof(fluid_particle_t*);
    total_bytes+=bytes;
    fluid_sim->fluid_particle_pointers = malloc(bytes);
    if(fluid_sim->fluid_particle_pointers == NULL)
        printf("Could not allocate fluid_particle_pointers\n");

    // Allocate neighbor array
    fluid_sim->neighbor_grid->neighbors = calloc(fluid_sim->params->max_fluid_particles_local, sizeof(neighbor_t));
    fluid_particle_t **fluid_neighbors = calloc(fluid_sim->params->max_fluid_particles_local * fluid_sim->neighbor_grid->max_neighbors, sizeof(fluid_particle_t *));
    // Set pointer in each bucket
    for(i=0; i< fluid_sim->params->max_fluid_particles_local; i++ )
        fluid_sim->neighbor_grid->neighbors[i].fluid_neighbors = &(fluid_neighbors[i*fluid_sim->neighbor_grid->max_neighbors]);
    total_bytes+= (fluid_sim->params->max_fluid_particles_local*sizeof(neighbor_t) + fluid_sim->neighbor_grid->max_neighbors*sizeof(fluid_particle_t *));
    if(fluid_sim->neighbor_grid->neighbors == NULL || fluid_neighbors == NULL)
        printf("Could not allocate neighbors\n");

    // UNIFORM GRID HASH
    fluid_sim->neighbor_grid->size_x = ceil((fluid_sim->boundary_global->max_x - fluid_sim->boundary_global->min_x) 
                                       / fluid_sim->neighbor_grid->spacing);
    fluid_sim->neighbor_grid->size_y = ceil((fluid_sim->boundary_global->max_y - fluid_sim->boundary_global->min_y) 
                                       / fluid_sim->neighbor_grid->spacing);
    unsigned int length_hash = fluid_sim->neighbor_grid->size_x * fluid_sim->neighbor_grid->size_y;
    printf("grid x: %d grid y %d\n", fluid_sim->neighbor_grid->size_x, fluid_sim->neighbor_grid->size_y);
    fluid_sim->neighbor_grid->grid_buckets = calloc(length_hash, sizeof(bucket_t));
    fluid_particle_t **bucket_particles = calloc(length_hash * fluid_sim->neighbor_grid->max_bucket_size, sizeof(fluid_particle_t *));
    for(i=0; i < length_hash; i++)
        fluid_sim->neighbor_grid->grid_buckets[i].fluid_particles = &(bucket_particles[i*fluid_sim->neighbor_grid->max_bucket_size]);
    total_bytes+= (length_hash * sizeof(bucket_t) + fluid_sim->neighbor_grid->max_bucket_size * sizeof(fluid_particle_t *));
    if(fluid_sim->neighbor_grid->grid_buckets == NULL || bucket_particles == NULL)
        printf("Could not allocate hash\n");

    // Allocate edge index arrays
    fluid_sim->edges->edge_pointers_left = malloc(fluid_sim->edges->max_edge_particles * sizeof(fluid_particle_t*));
    fluid_sim->edges->edge_pointers_right = malloc(fluid_sim->edges->max_edge_particles * sizeof(fluid_particle_t*));
    // Allocate out of bound index arrays
    fluid_sim->out_of_bounds->oob_pointer_indicies_left = malloc(fluid_sim->out_of_bounds->max_oob_particles * sizeof(int));
    fluid_sim->out_of_bounds->oob_pointer_indicies_right = malloc(fluid_sim->out_of_bounds->max_oob_particles * sizeof(int));
    fluid_sim->out_of_bounds->vacant_indicies = malloc(2*fluid_sim->out_of_bounds->max_oob_particles * sizeof(int));

    printf("bytes allocated: %lu\n", total_bytes);
}

// Free main simulation memory
void free_sim_memory(fluid_sim_t *fluid_sim)
{
    free(fluid_sim->fluid_particles);
    free(fluid_sim->fluid_particle_coords);
    free(fluid_sim->fluid_particle_pointers);
    free(fluid_sim->neighbor_grid->neighbors[0].fluid_neighbors);
    free(fluid_sim->neighbor_grid->neighbors);
    free(fluid_sim->neighbor_grid->grid_buckets[0].fluid_particles);
    free(fluid_sim->neighbor_grid->grid_buckets);
    free(fluid_sim->edges->edge_pointers_left);
    free(fluid_sim->edges->edge_pointers_right);
    free(fluid_sim->out_of_bounds->oob_pointer_indicies_left);
    free(fluid_sim->out_of_bounds->oob_pointer_indicies_right);
    free(fluid_sim->out_of_bounds->vacant_indicies);
}

// Allocate base structs used for simulation
void alloc_sim_structs(fluid_sim_t *fluid_sim)
{
   fluid_sim->params = (param_t*) calloc(1, sizeof(param_t));
   fluid_sim->water_volume_global = (AABB_t*) calloc(1, sizeof(AABB_t));
   fluid_sim->boundary_global = (AABB_t*) calloc(1, sizeof(AABB_t));
   fluid_sim->edges = (edge_t*) calloc(1, sizeof(edge_t));
   fluid_sim->out_of_bounds = (oob_t*) calloc(1, sizeof(oob_t)); 
}

// Free structs used for simulation
void free_sim_structs(fluid_sim_t *fluid_sim)
{
    free(fluid_sim->params);
    free(fluid_sim->water_volume_global);
    free(fluid_sim->boundary_global);
    free(fluid_sim->edges);
    free(fluid_sim->out_of_bounds);
}

// Initialize fluid parameters
// Additionally set world boudnary
void init_params(fluid_sim_t *fluid_sim)
{
    param_t *params = fluid_sim->params;

    params->tunable_params.kill_sim = false;
    params->tunable_params.active = true;
    params->tunable_params.g = 10.0f;
    params->tunable_params.time_step = 1.0f/60.0f;
    params->tunable_params.k = 0.2f;
    params->tunable_params.k_spring = 10.0f;
    params->tunable_params.sigma = 5.0f;
    params->tunable_params.beta = 0.5f;
    params->tunable_params.rest_density = 10.0;
    params->tunable_params.mover_width = 10.0f;
    params->tunable_params.mover_height = 10.0f;
    params->steps_per_frame = 4;  // Number of steps to compute before updating render node
    //params->tunable_params.time_step /= (float)steps_per_frame;

    // The number of particles used may differ slightly
    params->number_fluid_particles_global = 1500;

    // Boundary box
    // This simulation assumes in various spots min is 0.0
    fluid_sim->boundary_global->min_x = 0.0f;
    fluid_sim->boundary_global->max_x = 100.0f;
    fluid_sim->boundary_global->min_y = 0.0f;

    // Receive aspect ratio to scale world y max
    short pixel_dims[2];
    float aspect_ratio;
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
    aspect_ratio = (float)pixel_dims[0]/(float)pixel_dims[1];
    fluid_sim->boundary_global->max_y = fluid_sim->boundary_global->max_x / aspect_ratio;

    // water volume
    fluid_sim->water_volume_global->min_x = 10.0f;
    fluid_sim->water_volume_global->max_x = fluid_sim->boundary_global->max_x-10.0f;
    fluid_sim->water_volume_global->min_y = 10.0f;
    fluid_sim->water_volume_global->max_y = fluid_sim->boundary_global->max_y-10.0f;

    params->number_halo_particles = 0;
    params->number_halo_particles_left = 0;
    params->number_halo_particles_right = 0;
}

void construct_fluid_volume(fluid_sim_t *fluid_sim, float start_x, int number_particles_x)
{
    int num_y;

    // Unpack fluid_sim
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
    fluid_particle_t *fluid_particles = fluid_sim->fluid_particles;
    AABB_t* fluid = fluid_sim->water_volume_global;
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    float spacing = fluid_sim->params->tunable_params.smoothing_radius/2.0f;   

    // Number of particles in y,z, number in x is passed in
    num_y = floor((fluid->max_y - fluid->min_y ) / spacing);
    
    // zero out number of edge particles
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;

    // Place particles inside bounding volume
    float x,y;
    int nx = 0;
    int ny = 0;
    int i = 0;
    fluid_particle_t *p;

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

void partition_simulation(fluid_sim_t *fluid_sim, float *start_x, int *number_particles_x)
{
    // Fluid area in initial configuration
    float area = (fluid_sim->water_volume_global->max_x - fluid_sim->water_volume_global->min_x) 
                 * (fluid_sim->water_volume_global->max_y - fluid_sim->water_volume_global->min_y);

    // Initial spacing between particles
    float spacing_particle = pow(area/fluid_sim->params->number_fluid_particles_global,1.0/2.0);

    // Divide problem set amongst nodes
    partition_geometry(fluid_sim, start_x, number_particles_x, spacing_particle);

    // Set local/global number of particles to allocate
    set_particle_numbers(fluid_sim, *number_particles_x, spacing_particle);

    // We will allocate enough room for all particles on single node
    // We also must take into account halo particles are placed onto the end of the max particle index
    // So this value can be even greater than the number of global
    // Before reaching this point the program should, but doesn't, intelligenly clean up fluid_particles
    fluid_sim->params->max_fluid_particles_local = 2*fluid_sim->params->number_fluid_particles_global;

    // Smoothing radius, h
    fluid_sim->params->tunable_params.smoothing_radius = 2.0f*spacing_particle;

    // Particle mass is used to make density particle number inconsistant
    fluid_sim->params->particle_mass = (area*fluid_sim->params->tunable_params.rest_density)/(float)fluid_sim->params->number_fluid_particles_global;

    // Send initial world dimensions and max particle count to render node
    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    if(rank == 0) {
        float world_dims[2];
        world_dims[0] = fluid_sim->boundary_global->max_x;
        world_dims[1] = fluid_sim->boundary_global->max_y;
        MPI_Send(world_dims, 2, MPI_FLOAT, 0, 8, MPI_COMM_WORLD);
        MPI_Send(&fluid_sim->params->number_fluid_particles_global, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }
}

// Sets upper bound on number of particles, used for memory allocation
// These numbers are set judiciously for TitanTitan as the number of particles is always small
void set_particle_numbers(fluid_sim_t *fluid_sim, int number_particles_x, float spacing)
{
    int num_x, num_y, max_y;

    // Unpack fluid_sim
    AABB_t* boundary_global = fluid_sim->boundary_global;
    AABB_t* fluid_global = fluid_sim->water_volume_global;
    edge_t *edges = fluid_sim->edges;
    oob_t *out_of_bounds = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

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

    out_of_bounds->number_vacancies = 0;
}

// Set local boundary and fluid particle
void partition_geometry(fluid_sim_t *fluid_sim, float *x_start, int *num_particles_x, float spacing)
{

    // Unpack fluid_sim
    AABB_t *fluid_global = fluid_sim->water_volume_global;
    AABB_t *boundary_global = fluid_sim->boundary_global;
    param_t *params = fluid_sim->params;

    int i, rank, nprocs;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
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
    *num_particles_x = particle_length_x[rank];
        
    // Set node partition values
    params->tunable_params.node_start_x = fluid_global->min_x + ((number_to_left-1) * spacing);
    params->tunable_params.node_end_x   = params->tunable_params.node_start_x + (particle_length_x[rank] * spacing);
    
    if (rank == 0)
        params->tunable_params.node_start_x  = boundary_global->min_x;
    if (rank == nprocs-1)
        params->tunable_params.node_end_x   = boundary_global->max_x;

    printf("Rank %d start_x: %f, end_x :%f\n", rank, params->tunable_params.node_start_x, params->tunable_params.node_end_x);

    // Update requested number of particles with actual value used
    int num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    int total_x = 0;
    for(i=0; i<nprocs; i++)
        total_x += particle_length_x[i];
    params->number_fluid_particles_global = total_x * num_y;

    free(particle_length_x);
}

// initial Syncronization of tunable params with render node
void sync_initial_params(fluid_sim_t *fluid_sim)
{
    // Send intiial paramaters to render node
    tunable_parameters_t *null_tunable_param = NULL;
    int *null_recvcnts = NULL;
    int *null_displs = NULL;
    MPI_Gatherv(&fluid_sim->params->tunable_params, 1, TunableParamtype, null_tunable_param, null_recvcnts, null_displs, TunableParamtype, 0, MPI_COMM_WORLD);
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
