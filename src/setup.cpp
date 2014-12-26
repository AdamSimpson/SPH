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
#include "hash_sort.h"

// Main memory allocation and particle initialization
void alloc_and_init_sim(fluid_sim_t *fluid_sim)
{
    float start_x = 0.0f;
    int number_particles_x = 0;

    // Partition simulation and set particle numbers
    // Requires MPI_Send for particle counts
    partition_simulation(fluid_sim, &start_x, &number_particles_x);

    // Allocate main simulation memory and set struct values
    alloc_sim(fluid_sim);

    // Initialize particles
    init_sim_particles(fluid_sim, start_x, number_particles_x);

    // Print some parameters
    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, fluid_sim->params->number_fluid_particles_local, fluid_sim->params->tunable_params.smoothing_radius);
}

// Allocate main simulation memory
void alloc_sim(fluid_sim_t *fluid_sim)
{
    unsigned int i;

    // Neighbor grid setup
    fluid_sim->neighbor_grid = (neighbor_grid_t*)malloc(sizeof(neighbor_grid_t));
    if(fluid_sim->neighbor_grid == NULL)
        printf("Could not allocate neighbor grid\n");
    fluid_sim->neighbor_grid->max_neighbors = 50;
    fluid_sim->neighbor_grid->spacing = fluid_sim->params->tunable_params.smoothing_radius;

    size_t total_bytes = 0;
    size_t bytes;

    // Allocate fluid particles arrays
    bytes = fluid_sim->params->max_fluid_particles_local * sizeof(float);

    total_bytes+=bytes;
    fluid_sim->fluid_particles->x_star = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->x_star == NULL)
        printf("Could not allocate x_star\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->y_star = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->y_star == NULL)
        printf("Could not allocate y_star\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->z_star = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->z_star == NULL)
        printf("Could not allocate z_star\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->x = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->x == NULL)
        printf("Could not allocate x\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->y = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->y == NULL)
        printf("Could not allocate y\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->z = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->z == NULL)
        printf("Could not allocate z\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->v_x = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->v_x == NULL)
        printf("Could not allocate v_x\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->v_y = (float*) malloc(bytes);

    if(fluid_sim->fluid_particles->v_y == NULL)
        printf("Could not allocate v_y\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->v_z = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->v_z == NULL)
        printf("Could not allocate v_z\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->dp_x = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->dp_x == NULL)
        printf("Could not allocate dp_x\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->dp_y = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->dp_y == NULL)
        printf("Could not allocate dp_y\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->dp_z = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->dp_z == NULL)
        printf("Could not allocate dp_z\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->density = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->density == NULL)
        printf("Could not allocate density\n");
    total_bytes+=bytes;
    fluid_sim->fluid_particles->lambda = (float*) malloc(bytes);
    if(fluid_sim->fluid_particles->lambda == NULL)
        printf("Could not allocate lambda\n");
    bytes = fluid_sim->params->max_fluid_particles_local * sizeof(uint);
    total_bytes+=bytes;
    fluid_sim->fluid_particles->id = (uint*) malloc(bytes);
    if(fluid_sim->fluid_particles->id == NULL) 
        printf("Could not allocate id\n");


    // Allocate (x,y,z) coordinate array, transfer pixel coords
    bytes = 3 * fluid_sim->params->max_fluid_particles_local * sizeof(short);
    total_bytes+=bytes;
    fluid_sim->fluid_particle_coords = (short*)malloc(bytes);
    if(fluid_sim->fluid_particle_coords == NULL)
        printf("Could not allocate fluid_particle coords\n");

    // Allocate index array used to traverse non vacant particles
    bytes = fluid_sim->params->max_fluid_particles_local * sizeof(uint);
    total_bytes+=bytes;
    fluid_sim->fluid_particle_indices = (uint*)malloc(bytes);
    if(fluid_sim->fluid_particle_indices == NULL)
        printf("Could not allocate fluid_particle_indices\n");

    // Allocate neighbor array
    fluid_sim->neighbor_grid->neighbors = (neighbor_t*)calloc(fluid_sim->params->max_fluid_particles_local, sizeof(neighbor_t));
    uint *fluid_neighbors = (uint*)calloc(fluid_sim->params->max_fluid_particles_local * fluid_sim->neighbor_grid->max_neighbors, sizeof(uint));
    // Set pointer in each bucket
    for(i=0; i< fluid_sim->params->max_fluid_particles_local; i++ )
        fluid_sim->neighbor_grid->neighbors[i].fluid_neighbors = &(fluid_neighbors[i*fluid_sim->neighbor_grid->max_neighbors]);
    total_bytes+= (fluid_sim->params->max_fluid_particles_local*sizeof(neighbor_t) + fluid_sim->neighbor_grid->max_neighbors*sizeof(uint));
    if(fluid_sim->neighbor_grid->neighbors == NULL || fluid_neighbors == NULL)
        printf("Could not allocate neighbors\n");

    printf("max %f, min %f\n", fluid_sim->boundary_global->max_x, fluid_sim->boundary_global->min_x);


    // UNIFORM GRID HASH
    fluid_sim->neighbor_grid->size_x = ceil((fluid_sim->boundary_global->max_x - fluid_sim->boundary_global->min_x) 
                                       / fluid_sim->neighbor_grid->spacing);
    fluid_sim->neighbor_grid->size_y = ceil((fluid_sim->boundary_global->max_y - fluid_sim->boundary_global->min_y) 
                                       / fluid_sim->neighbor_grid->spacing);
    fluid_sim->neighbor_grid->size_z = ceil((fluid_sim->boundary_global->max_z - fluid_sim->boundary_global->min_z)
                                       / fluid_sim->neighbor_grid->spacing);
    unsigned int length_hash = fluid_sim->neighbor_grid->size_x * fluid_sim->neighbor_grid->size_y *  fluid_sim->neighbor_grid->size_z;
    printf("grid x: %d grid y: %d grid z: %d\n", fluid_sim->neighbor_grid->size_x, fluid_sim->neighbor_grid->size_y,  fluid_sim->neighbor_grid->size_z);

    // Start index for hash values
    fluid_sim->neighbor_grid->start_indices = (uint*)calloc(length_hash, sizeof(uint));

    // End index for hash values
    fluid_sim->neighbor_grid->end_indices = (uint*)calloc(length_hash, sizeof(uint));
    // Array of hash values
    fluid_sim->neighbor_grid->hash_values = (uint*)calloc(fluid_sim->params->max_fluid_particles_local, sizeof(uint));

    // Array of particle id's
    fluid_sim->neighbor_grid->particle_ids = (uint*)calloc(fluid_sim->params->max_fluid_particles_local, sizeof(uint));

    total_bytes+= (length_hash + fluid_sim->params->max_fluid_particles_local) * sizeof(uint);
    if(fluid_sim->neighbor_grid->start_indices  == NULL || fluid_sim->neighbor_grid->end_indices == NULL)
        printf("Could not allocate start/end indices\n");
    if(fluid_sim->neighbor_grid->hash_values  == NULL || fluid_sim->neighbor_grid->particle_ids == NULL)
        printf("Could not allocate hash_values/particle_ids\n");

    // Allocate edge index arrays
    fluid_sim->edges->edge_indices_left = (uint*)malloc(fluid_sim->edges->max_edge_particles * sizeof(uint));
    fluid_sim->edges->edge_indices_right = (uint*)malloc(fluid_sim->edges->max_edge_particles * sizeof(uint));
    // Allocate out of bound index arrays
    fluid_sim->out_of_bounds->oob_index_indices_left = (uint*)malloc(fluid_sim->out_of_bounds->max_oob_particles * sizeof(uint));
    fluid_sim->out_of_bounds->oob_index_indices_right = (uint*)malloc(fluid_sim->out_of_bounds->max_oob_particles * sizeof(uint));
    fluid_sim->out_of_bounds->vacant_indices = (uint*)malloc(2*fluid_sim->out_of_bounds->max_oob_particles * sizeof(uint));

    printf("bytes allocated: %lu\n", total_bytes);
}

// Free main simulation memory
void free_sim_memory(fluid_sim_t *fluid_sim)
{
    free(fluid_sim->fluid_particles->x_star);
    free(fluid_sim->fluid_particles->y_star);
    free(fluid_sim->fluid_particles->z_star);
    free(fluid_sim->fluid_particles->x);
    free(fluid_sim->fluid_particles->y);
    free(fluid_sim->fluid_particles->z);
    free(fluid_sim->fluid_particles->v_x);
    free(fluid_sim->fluid_particles->v_y);
    free(fluid_sim->fluid_particles->v_z);
    free(fluid_sim->fluid_particles->dp_x);
    free(fluid_sim->fluid_particles->dp_y);
    free(fluid_sim->fluid_particles->dp_z);
    free(fluid_sim->fluid_particles->density);
    free(fluid_sim->fluid_particles->lambda);
    free(fluid_sim->fluid_particles->id);
    free(fluid_sim->fluid_particle_coords);
    free(fluid_sim->fluid_particle_indices);
    free(fluid_sim->neighbor_grid->neighbors[0].fluid_neighbors);
    free(fluid_sim->neighbor_grid->neighbors);
    free(fluid_sim->edges->edge_indices_left);
    free(fluid_sim->edges->edge_indices_right);
    free(fluid_sim->out_of_bounds->oob_index_indices_left);
    free(fluid_sim->out_of_bounds->oob_index_indices_right);
    free(fluid_sim->out_of_bounds->vacant_indices);
}

void init_sim_particles(fluid_sim_t *fluid_sim, float start_x, int number_particles_x)
{
    int i;
    uint p_index;

    // Create fluid volume
    construct_fluid_volume(fluid_sim, start_x, number_particles_x);

    // invalidate unused fluid pointers
    for(i=fluid_sim->params->number_fluid_particles_local; i<fluid_sim->params->max_fluid_particles_local; i++)
        fluid_sim->fluid_particle_indices[i] = ((uint)-1);

    // Initialize particle values
    for(i=0; i<fluid_sim->params->number_fluid_particles_local; i++) {
        fluid_sim->fluid_particles->v_x[i] = 0.0f;
        fluid_sim->fluid_particles->v_y[i] = 0.0f;
        fluid_sim->fluid_particles->v_z[i] = 0.0f;
    }
}

// Allocate base structs used for simulation
void alloc_sim_structs(fluid_sim_t *fluid_sim)
{
   fluid_sim->params = (param_t*) calloc(1, sizeof(param_t));
   fluid_sim->fluid_particles = (fluid_particles_t*) calloc(1, sizeof(fluid_particles_t));
   fluid_sim->water_volume_global = (AABB_t*) calloc(1, sizeof(AABB_t));
   fluid_sim->boundary_global = (AABB_t*) calloc(1, sizeof(AABB_t));
   fluid_sim->edges = (edge_t*) calloc(1, sizeof(edge_t));
   fluid_sim->out_of_bounds = (oob_t*) calloc(1, sizeof(oob_t)); 
}

// Free structs used for simulation
void free_sim_structs(fluid_sim_t *fluid_sim)
{
    free(fluid_sim->params);
    free(fluid_sim->fluid_particles);
    free(fluid_sim->water_volume_global);
    free(fluid_sim->boundary_global);
    free(fluid_sim->edges);
    free(fluid_sim->out_of_bounds);
}

// Initialize fluid parameters
// Additionally set world boudnary
void init_params(fluid_sim_t *fluid_sim)
{
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    params->tunable_params.kill_sim = false;
    params->tunable_params.active = true;
    params->tunable_params.g = 10.0f;
    params->tunable_params.time_step = 1.0f/60.0f;
    params->tunable_params.k = 0.2f;
    params->tunable_params.c = 0.01;
    params->tunable_params.rest_density = 0.1;
    params->tunable_params.mover_radius = 25.0f;
    params->steps_per_frame = 4;  // Number of steps to compute before updating render node
    //params->tunable_params.time_step /= (float)steps_per_frame;

    // The number of particles used may differ slightly
    params->number_fluid_particles_global = 5000;

    // Boundary box
    // This simulation assumes in various spots min is 0.0
    fluid_sim->boundary_global->min_x = 0.0f;
    fluid_sim->boundary_global->max_x = 100.0f;
    fluid_sim->boundary_global->min_y = 0.0f;
    fluid_sim->boundary_global->max_y = 40.0f;
    fluid_sim->boundary_global->min_z = 0.0f;
    fluid_sim->boundary_global->max_z = 30.0f;

    // Receive aspect ratio to scale world y max
    short pixel_dims[2];
    float aspect_ratio;
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
    aspect_ratio = (float)pixel_dims[0]/(float)pixel_dims[1];
//    fluid_sim->boundary_global->max_z = fluid_sim->boundary_global->max_x / aspect_ratio;

    // water volume
    fluid_sim->water_volume_global->min_x = 10.0f;
    fluid_sim->water_volume_global->max_x = fluid_sim->boundary_global->max_x-10.0f;
    fluid_sim->water_volume_global->min_y = 15.0f;
    fluid_sim->water_volume_global->max_y = fluid_sim->boundary_global->max_y-15.0f;
    fluid_sim->water_volume_global->min_z = 5.0f;
    fluid_sim->water_volume_global->max_z = fluid_sim->boundary_global->max_z-5.0f;

    // Zero out number of halo particles
    params->number_halo_particles = 0;
    params->number_halo_particles_left = 0;
    params->number_halo_particles_right = 0;

    // zero out number of edge particles
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;
}

void construct_fluid_volume(fluid_sim_t *fluid_sim, float start_x, int number_particles_x)
{
    int num_y,num_z;

    // Unpack fluid_sim
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;

    AABB_t* fluid = fluid_sim->water_volume_global;
    param_t *params = fluid_sim->params;

    float spacing = fluid_sim->params->tunable_params.smoothing_radius/2.0f;   

    // Number of particles in y,z, number in x is passed in
    num_y = floor((fluid->max_y - fluid->min_y ) / spacing);
    num_z = floor((fluid->max_z - fluid->min_z ) / spacing);    

    // Place particles inside bounding volume
    float x,y,z;
    int nx = 0;
    int ny = 0;
    int nz = 0;
    int i = 0;
    uint p;

    for(nz=0; nz<num_z; nz++) {
            z = fluid->min_z + nz*spacing;
        for(ny=0; ny<num_y; ny++) {
            y = fluid->min_y + ny*spacing;
            for(nx=0; nx<number_particles_x; nx++) {
                x = fluid->min_x + (start_x + nx)*spacing;
                fluid_particles->x[i] = x;
                fluid_particles->y[i] = y;
                fluid_particles->z[i] = z;            

                // Set index array
                fluid_particle_indices[i] = i;
	        fluid_particles->id[i] = i;
                i++;
            }
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
    // Fluid volume in initial configuration
    float volume = (fluid_sim->water_volume_global->max_x - fluid_sim->water_volume_global->min_x) 
                 * (fluid_sim->water_volume_global->max_y - fluid_sim->water_volume_global->min_y)
                 * (fluid_sim->water_volume_global->max_z - fluid_sim->water_volume_global->min_z);

    // Initial spacing between particles
    float spacing_particle = pow(volume/fluid_sim->params->number_fluid_particles_global,1.0/3.0);

    // Divide problem set amongst nodes
    partition_geometry(fluid_sim, start_x, number_particles_x, spacing_particle);

    // Set local/global number of particles to allocate
    set_particle_numbers(fluid_sim, *number_particles_x, spacing_particle);

    // We also must take into account halo particles are placed onto the end of the max particle index
    // So this value can be even greater than the number of global
    // Before reaching this point the program should, but doesn't, intelligenly clean up fluid_particles
    fluid_sim->params->max_fluid_particles_local = 2*fluid_sim->params->number_fluid_particles_global;

    // Smoothing radius, h
    fluid_sim->params->tunable_params.smoothing_radius = 2.0f*spacing_particle;
    // Set dq, used in S_corr
    fluid_sim->params->tunable_params.dq = 0.3*fluid_sim->params->tunable_params.smoothing_radius;


    // Particle mass is used to make density particle number inconsistant
    fluid_sim->params->particle_mass = (volume*fluid_sim->params->tunable_params.rest_density)/(float)fluid_sim->params->number_fluid_particles_global;

    // Send initial world dimensions and max particle count to render node
    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    if(rank == 0) {
        float world_dims[3];
        world_dims[0] = fluid_sim->boundary_global->max_x;
        world_dims[1] = fluid_sim->boundary_global->max_y;
        world_dims[2] = fluid_sim->boundary_global->max_z;
        MPI_Send(world_dims, 3, MPI_FLOAT, 0, 8, MPI_COMM_WORLD);
        MPI_Send(&fluid_sim->params->number_fluid_particles_global, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }

}

// Sets upper bound on number of particles, used for memory allocation
// These numbers are set judiciously for TitanTitan as the number of particles is always small
void set_particle_numbers(fluid_sim_t *fluid_sim, int number_particles_x, float spacing)
{
    int num_x, num_y, num_z;

    // Unpack fluid_sim
    AABB_t* fluid_global = fluid_sim->water_volume_global;
    edge_t *edges = fluid_sim->edges;
    oob_t *out_of_bounds = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

    // Set fluid local
    num_x = number_particles_x;
    num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    num_z = floor((fluid_global->max_z - fluid_global->min_z ) / spacing);

    // Maximum edge(halo) particles
    edges->max_edge_particles = params->number_fluid_particles_global;

    // The out of bounds particles can become quite large
    // If a flood of particles flows into and then out of a node
    // This will be large
    out_of_bounds->max_oob_particles = params->number_fluid_particles_global;

    // Initial fluid particles
    int num_initial = num_x * num_y * num_z;
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
    int *particle_length_x = (int*)malloc(nprocs*sizeof(int));
    
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
    params->tunable_params.proc_start = fluid_global->min_x + ((number_to_left-1) * spacing);
    params->tunable_params.proc_end   = params->tunable_params.proc_start + (particle_length_x[rank] * spacing);
    
    if (rank == 0)
        params->tunable_params.proc_start  = boundary_global->min_x;
    if (rank == nprocs-1)
        params->tunable_params.proc_end   = boundary_global->max_x;

    printf("Rank %d start_x: %f, end_x :%f\n", rank, params->tunable_params.proc_start, params->tunable_params.proc_end);

    // Update requested number of particles with actual value used
    int num_y = floor((fluid_global->max_y - fluid_global->min_y ) / spacing);
    int num_z = floor((fluid_global->max_z - fluid_global->min_z ) / spacing);
    int total_x = 0;
    for(i=0; i<nprocs; i++)
        total_x += particle_length_x[i];

    params->number_fluid_particles_global = total_x * num_y * num_z;

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
