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
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "mpi.h"
#include "hash.h"
#include "renderer.h"
#include "geometry.h"
#include "fluid.h"
#include "communication.h"

#include "cuda.h"

#ifdef LIGHT
#include "rgb_light.h"
#include <unistd.h>
#endif

int main(int argc, char *argv[])
{
    int return_value;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank;

    // Rank in world space
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    create_communicators();

    createMpiTypes();

    // Rank 0 is the render node, otherwise a simulation node
    if(rank == 0)
        return_value = start_renderer();
    else
        start_simulation();

    MPI_Finalize();
    return return_value;
}

void start_simulation()
{
    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    printf("compute rank: %d, num compute procs: %d \n",rank, nprocs);

    param *params;
    cudaMallocManaged(&params, sizeof(param));
    AABB_t water_volume_global;
    AABB_t *boundary_global;
    cudaMallocManaged(&boundary_global, sizeof(AABB_t));
    edge_t edges;
    oob_t out_of_bounds;

    unsigned int i;

    params->tunable_params.kill_sim = false;
    params->tunable_params.active = true;
    params->tunable_params.g = 6.0f;
    params->tunable_params.time_step = 1.0f/35.0f/4.0f;
    params->tunable_params.k = 0.2f;
    params->tunable_params.k_near = 6.0f;
    params->tunable_params.k_spring = 10.0f;
    params->tunable_params.sigma = 5.0f;
    params->tunable_params.beta = 0.5f;
    params->tunable_params.rest_density = 30.0f;
    params->tunable_params.mover_width = 2.0f;
    params->tunable_params.mover_height = 2.0f;
    params->tunable_params.mover_type = SPHERE_MOVER;

    // The number of particles used may differ slightly
    #ifdef RASPI
    params->number_fluid_particles_global = 1500;
    #else
    params->number_fluid_particles_global = 5500;
    #endif

    // Boundary box
    // This simulation assumes in various spots min is 0.0
    boundary_global->min_x = 0.0f;
    boundary_global->max_x = 15.0f;
    boundary_global->min_y = 0.0f;

    // Receive aspect ratio to scale world y max
    short pixel_dims[2];
    float aspect_ratio;
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
    aspect_ratio = (float)pixel_dims[0]/(float)pixel_dims[1];
    boundary_global->max_y = boundary_global->max_x / aspect_ratio;

    // water volume
    water_volume_global.min_x = 0.0f;
    water_volume_global.max_x = boundary_global->max_x;
    water_volume_global.min_y = 0.0f;
    water_volume_global.max_y = boundary_global->max_y;

    params->number_halo_particles = 0;

    int start_x;  // where in x direction this nodes particles start
    int number_particles_x; // number of particles in x direction for this node

    // Fluid area in initial configuration
    float area = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y);

    // Initial spacing between particles
    float spacing_particle = pow(area/params->number_fluid_particles_global,1.0/2.0);

    // Divide problem set amongst nodes
    partitionProblem(boundary_global, &water_volume_global, &start_x, &number_particles_x, spacing_particle, params);

    // Set local/global number of particles to allocate
    setParticleNumbers(boundary_global, &water_volume_global, &edges, &out_of_bounds, number_particles_x, spacing_particle, params);

    // We will allocate enough room for all particles on single node
    // We also must take into account halo particles are placed onto the end of the max particle index
    // So this value can be even greater than the number of global
    // Before reaching this point the program should, but doesn't, intelligenly clean up fluid_particles
    int max_fluid_particles_local = 2*params->number_fluid_particles_global;

    // Smoothing radius, h
    params->tunable_params.smoothing_radius = 2.0f*spacing_particle;

    printf("smoothing radius: %f\n", params->tunable_params.smoothing_radius);

    // Send initial world dimensions and max particle count to render node
    if(rank == 0) {
        float world_dims[2];
        world_dims[0] = boundary_global->max_x;
        world_dims[1] = boundary_global->max_y;
        MPI_Send(world_dims, 2, MPI_FLOAT, 0, 8, MPI_COMM_WORLD);
	MPI_Send(&params->number_fluid_particles_global, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }

    // Neighbor grid setup
//    neighbor_grid_t neighbor_grid;
//    neighbor_grid.max_bucket_size = 100;
//    neighbor_grid.max_neighbors = neighbor_grid.max_bucket_size*4;
    params->grid_spacing = params->tunable_params.smoothing_radius;

    size_t total_bytes = 0;
    size_t bytes;
    // Allocate fluid particles array
    bytes = max_fluid_particles_local * sizeof(fluid_particle);
    total_bytes+=bytes;
    fluid_particle *fluid_particles = malloc(bytes);
    if(fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");

    // Allocate (x,y) coordinate array, transfer pixel coords
    bytes = 2 * max_fluid_particles_local * sizeof(short);
    total_bytes+=bytes;
    short *fluid_particle_coords = malloc(bytes);
    if(fluid_particle_coords == NULL)
        printf("Could not allocate fluid_particle coords\n");

    // Allocate pointer array used to traverse non vacant particles
    bytes = max_fluid_particles_local * sizeof(fluid_particle*);
    total_bytes+=bytes;
    fluid_particle **fluid_particle_pointers;
    cudaMallocManaged(&fluid_particle_pointers ,bytes);
    if(fluid_particle_pointers == NULL)
        printf("Could not allocate fluid_particle_pointers\n");

    // Allocate uint hash value array
    bytes = max_fluid_particles_local * sizeof(unsigned int);
    total_bytes += bytes;
    unsigned int * hash_values;
    cudaMalloc(&hash_values, bytes);    

    // Allocate uint particle ID array
    bytes = max_fluid_particles_local * sizeof(unsigned int);
    total_bytes += bytes;
    unsigned int *particle_ids;
    cudaMalloc(&particle_ids, bytes);

    params->grid_size_x = ceil((boundary_global->max_x - boundary_global->min_x) / params->grid_spacing);
    params->grid_size_y = ceil((boundary_global->max_y - boundary_global->min_y) / params->grid_spacing);
    unsigned int length_hash = params->grid_size_x * params->grid_size_y;

    // Allocate uint cell start index array
    bytes = length_hash * sizeof(unsigned int);
    total_bytes += bytes;
    unsigned int * start_indexes;
    cudaMalloc(&start_indexes, bytes);

    // Allocate uint cell end index array
    bytes = length_hash * sizeof(unsigned int);
    total_bytes += bytes;
    unsigned int *end_indexes;
    cudaMalloc(&end_indexes, bytes);

/*
    // Allocate neighbor array
    neighbor *neighbors = calloc(max_fluid_particles_local, sizeof(neighbor));
    fluid_particle **fluid_neighbors = calloc(max_fluid_particles_local * neighbor_grid.max_neighbors, sizeof(fluid_particle *));
    // Set pointer in each bucket
    for(i=0; i< max_fluid_particles_local; i++ )
        neighbors[i].fluid_neighbors = &(fluid_neighbors[i*neighbor_grid.max_neighbors]);

    neighbor_grid.neighbors = neighbors;
    total_bytes+= (max_fluid_particles_local*sizeof(neighbor) + neighbor_grid.max_neighbors*sizeof(fluid_particle *));
    if(neighbors == NULL || fluid_neighbors == NULL)
        printf("Could not allocate neighbors\n");
*/

/*
    // UNIFORM GRID HASH
    neighbor_grid.size_x = ceil((boundary_global.max_x - boundary_global.min_x) / neighbor_grid.spacing);
    neighbor_grid.size_y = ceil((boundary_global.max_y - boundary_global.min_y) / neighbor_grid.spacing);
    unsigned int length_hash = neighbor_grid.size_x * neighbor_grid.size_y;
    printf("grid x: %d grid y %d\n", neighbor_grid.size_x, neighbor_grid.size_y);
    bucket_t* grid_buckets = calloc(length_hash, sizeof(bucket_t));
    fluid_particle **bucket_particles = calloc(length_hash * neighbor_grid.max_bucket_size, sizeof(fluid_particle *));
    neighbor_grid.grid_buckets = grid_buckets;
    for(i=0; i < length_hash; i++)
	grid_buckets[i].fluid_particles = &(bucket_particles[i*neighbor_grid.max_bucket_size]);
    total_bytes+= (length_hash * sizeof(bucket_t) + neighbor_grid.max_bucket_size * sizeof(fluid_particle *));
    if(grid_buckets == NULL || bucket_particles == NULL)
        printf("Could not allocate hash\n");
*/

    // Allocate edge index arrays
    // Edge and out of bounds arrays are handled on host
    edges.edge_pointers_left = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    edges.edge_pointers_right = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    // Allocate out of bound index arrays
    out_of_bounds.oob_pointer_indicies_left = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.oob_pointer_indicies_right = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.vacant_indicies = malloc(2*out_of_bounds.max_oob_particles * sizeof(int));

    printf("bytes allocated: %lu\n", total_bytes);

    // Initialize particles
    initParticles(fluid_particle_pointers, fluid_particles, &water_volume_global, start_x,
		  number_particles_x, &edges, max_fluid_particles_local, spacing_particle, params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, params->number_fluid_particles_local, params->tunable_params.smoothing_radius);

    // Send intiial paramaters to render node
    tunable_parameters *null_tunable_param = NULL;
    int *null_recvcnts = NULL;
    int *null_displs = NULL;
    MPI_Gatherv(&params->tunable_params, 1, TunableParamtype, null_tunable_param, null_recvcnts, null_displs, TunableParamtype, 0, MPI_COMM_WORLD);

    // Initialize RGB Light if present
    #ifdef LIGHT
    rgb_light_t light_state;
    float *colors_by_rank = malloc(3*nprocs*sizeof(float));
    MPI_Bcast(colors_by_rank, 3*nprocs, MPI_FLOAT, 0, MPI_COMM_WORLD);
    init_rgb_light(&light_state, 255*colors_by_rank[3*rank], 255*colors_by_rank[3*rank+1], 255*colors_by_rank[3*rank+2]);
    free(colors_by_rank);
    // Without this pause the lights can sometimes change color too quickly the first time step
    sleep(1);
    #endif    

    fluid_particle *p;
    fluid_particle *null_particle = NULL;
    float *null_float = NULL;

    MPI_Request coords_req = MPI_REQUEST_NULL;

    #ifdef RASPI
    int steps_per_frame = 1; // Number of steps to compute before updating render node
    #else
    int steps_per_frame = 4; 
    #endif
    int sub_step = 0; // substep range from 0 to < steps_per_frame

    int threads_per_block = 256;
    int num_blocks;

    // Main simulation loop
    while(1) {

        // Initialize velocities
        apply_gravity(fluid_particle_pointers, params);

        // Viscosity impluse
        viscosity_impluses(fluid_particle_pointers, particle_ids, start_indexes, end_indexes, params);

        // Advance to predicted position and set OOB particles
        predict_positions(fluid_particle_pointers, boundary_global, params)

        // Make sure that async send to render node is complete
        if(sub_step == 0)
        {
            if(coords_req != MPI_REQUEST_NULL)
	        MPI_Wait(&coords_req, MPI_STATUS_IGNORE);
        }

        #ifdef LIGHT
        char previously_active = params->tunable_params.active;
        #endif

        // Receive updated paramaters from render nodes
        if(sub_step == steps_per_frame-1)
            MPI_Scatterv(null_tunable_param, 0, null_displs, TunableParamtype, params->tunable_params, 1, TunableParamtype, 0,  MPI_COMM_WORLD);

        #ifdef LIGHT
        // If recently added to computation turn light to light state color
        // If recently taken out of computation turn light to white
        char currently_active = params->tunable_params.active;
        if (!previously_active && currently_active)
            rgb_light_reset(&light_state);
        else if (!currently_active && previously_active)
            rgb_light_white(&light_state);
        #endif

        if(params->tunable_params.kill_sim)
            break;

        // Synchronize kernels
        cudaDeviceSynchrnoize();

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(fluid_particle_pointers, fluid_particles, &out_of_bounds, boundary_global, params);

         // Exchange halo particles
        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, params);
        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, params);

        hash_particles(fluid_particle_pointers, hash_values, particle_ids, start_indexes, end_indexes, params);

        // TO DO: update density kernel
        calculate_densities(fluid_particle_pointers, params);
        calculate_pressures(fluid_particle_pointers, params);

        // double density relaxation
        // halo particles will be missing origin contributions to density/pressure
        double_density_relaxation(fluid_particle_pointers, neighbors, params);

        // update velocity
        updateVelocities(fluid_particle_pointers, boundary_global, params);

        // Not updating halo particles and hash after relax can be used to speed things up
        // Not updating these can cause unstable behavior

        // Hash particles
        hash_particles(fluid_particle_pointers, hash_values, particle_ids, start_indexes, end_indexes, params);

        // Exchange halo particles from relaxed positions
        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, params);
        // Finish asynch halo exchange
        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, params);

        // We do not transfer particles that have gone OOB since relaxation
        // to reduce communication cost

        // Pack fluid particle coordinates
        // This sends results as short in pixel coordinates
        if(sub_step == steps_per_frame-1)
        {
            for(i=0; i<params->number_fluid_particles_local; i++) {
                p = fluid_particle_pointers[i];
                fluid_particle_coords[i*2] = (2.0f*p->x/boundary_global->max_x - 1.0f) * SHRT_MAX; // convert to short using full range
                fluid_particle_coords[(i*2)+1] = (2.0f*p->y/boundary_global->max_y - 1.0f) * SHRT_MAX; // convert to short using full range
            }
            // Async send fluid particle coordinates to render node
            MPI_Isend(fluid_particle_coords, 2*params->number_fluid_particles_local, MPI_SHORT, 0, 17, MPI_COMM_WORLD, &coords_req);
        }

        if(sub_step == steps_per_frame-1)
            sub_step = 0;
        else
	    sub_step++;

        // Synchronize kernels
        cudaDeviceSynchrnoize();

    }

    #ifdef LIGHT
    rgb_light_off(&light_state);
    #endif

    // Release memory
    free(fluid_particles);
    free(fluid_particle_coords);
    free(fluid_particle_pointers);
    free(neighbors);
    free(fluid_neighbors);
    free(grid_buckets);
    free(bucket_particles);
    free(edges.edge_pointers_left);
    free(edges.edge_pointers_right);
    free(out_of_bounds.oob_pointer_indicies_left);
    free(out_of_bounds.oob_pointer_indicies_right);
    free(out_of_bounds.vacant_indicies);

    // Close MPI
    freeMpiTypes();

}

// Identify out of bounds particles and send them to appropriate rank
void identify_oob_particles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, oob_t *out_of_bounds, AABB_t *boundary_global, param *params)
{
    int i;
    fluid_particle *p;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];

        // Set OOB particle indicies and update number
        if (p->x < params->tunable_params.node_start_x)
            out_of_bounds->oob_pointer_indicies_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (p->x > params->tunable_params.node_end_x)
            out_of_bounds->oob_pointer_indicies_right[out_of_bounds->number_oob_particles_right++] = i;
    }
 
   // Transfer particles that have left the processor bounds
   transferOOBParticles(fluid_particle_pointers, fluid_particles, out_of_bounds, params);
}

// Initialize particles
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                   AABB_t *water, int start_x, int number_particles_x,    
                   edge_t *edges, int max_fluid_particles_local, float spacing, param* params)
{
    int i;
    fluid_particle *p;

    // Create fluid volume
    constructFluidVolume(fluid_particle_pointers, fluid_particles, water, start_x, number_particles_x, edges, spacing, params);

    // NULL out unused fluid pointers
    for(i=params->number_fluid_particles_local; i<max_fluid_particles_local; i++)
        fluid_particle_pointers[i] = NULL;

    // Initialize particle values
    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particle_pointers[i]->a_x = 0.0f;
        fluid_particle_pointers[i]->a_y = 0.0f;
        fluid_particle_pointers[i]->v_x = 0.0f;
        fluid_particle_pointers[i]->v_y = 0.0f;
    }
}
