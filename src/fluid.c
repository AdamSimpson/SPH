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

#ifdef BLINK1
#include "blink1_light.h"
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

    param params;
    AABB_t water_volume_global;
    AABB_t boundary_global;
    edge_t edges;
    oob_t out_of_bounds;

    unsigned int i;

    params.tunable_params.kill_sim = false;
    params.tunable_params.active = true;
    params.tunable_params.g = 10.0f;
    params.tunable_params.time_step = 1.0f/60.0f;
    params.tunable_params.k = 0.2f;
    params.tunable_params.k_spring = 10.0f;
    params.tunable_params.sigma = 5.0f;
    params.tunable_params.beta = 0.5f;
    params.tunable_params.rest_density = 10.0;
    params.tunable_params.mover_width = 10.0f;
    params.tunable_params.mover_height = 10.0f;

    #ifdef RASPI
    int steps_per_frame = 1; // Number of steps to compute before updating render node
    #else
    int steps_per_frame = 4;
//    params.tunable_params.time_step /= (float)steps_per_frame;
    #endif

    // The number of particles used may differ slightly
    #ifdef RASPI
    params.number_fluid_particles_global = 1500;
    #else
    params.number_fluid_particles_global = 1500;
    #endif

    // Boundary box
    // This simulation assumes in various spots min is 0.0
    boundary_global.min_x = 0.0f;
    boundary_global.max_x = 100.0f;
    boundary_global.min_y = 0.0f;

    // Receive aspect ratio to scale world y max
    short pixel_dims[2];
    float aspect_ratio;
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
    aspect_ratio = (float)pixel_dims[0]/(float)pixel_dims[1];
    boundary_global.max_y = boundary_global.max_x / aspect_ratio;

    // water volume
    water_volume_global.min_x = 10.0f;
    water_volume_global.max_x = boundary_global.max_x-10.0f;
    water_volume_global.min_y = 10.0f;
    water_volume_global.max_y = boundary_global.max_y-10.0f;

    params.number_halo_particles = 0;
    params.number_halo_particles_left = 0;
    params.number_halo_particles_right = 0;

    int start_x;  // where in x direction this nodes particles start
    int number_particles_x; // number of particles in x direction for this node

    // Fluid area in initial configuration
    float area = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y);

    // Initial spacing between particles
    float spacing_particle = pow(area/params.number_fluid_particles_global,1.0/2.0);

    // Divide problem set amongst nodes
    partitionProblem(&boundary_global, &water_volume_global, &start_x, &number_particles_x, spacing_particle, &params);

    // Set local/global number of particles to allocate
    setParticleNumbers(&boundary_global, &water_volume_global, &edges, &out_of_bounds, number_particles_x, spacing_particle, &params);

    // We will allocate enough room for all particles on single node
    // We also must take into account halo particles are placed onto the end of the max particle index
    // So this value can be even greater than the number of global
    // Before reaching this point the program should, but doesn't, intelligenly clean up fluid_particles
    int max_fluid_particles_local = 2*params.number_fluid_particles_global;

    // Smoothing radius, h
    params.tunable_params.smoothing_radius = 2.0f*spacing_particle;

    printf("smoothing radius: %f\n", params.tunable_params.smoothing_radius);

    // Particle mass is used to make density particle number inconsistant
    params.particle_mass = (area*params.tunable_params.rest_density)/(float)params.number_fluid_particles_global;

    // Send initial world dimensions and max particle count to render node
    if(rank == 0) {
        float world_dims[2];
        world_dims[0] = boundary_global.max_x;
        world_dims[1] = boundary_global.max_y;
        MPI_Send(world_dims, 2, MPI_FLOAT, 0, 8, MPI_COMM_WORLD);
	MPI_Send(&params.number_fluid_particles_global, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }

    // Neighbor grid setup
    neighbor_grid_t neighbor_grid;
    neighbor_grid.max_bucket_size = 100;
    neighbor_grid.max_neighbors = neighbor_grid.max_bucket_size*4;
    neighbor_grid.spacing = params.tunable_params.smoothing_radius;

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
    fluid_particle **fluid_particle_pointers = malloc(bytes);
    if(fluid_particle_pointers == NULL)
        printf("Could not allocate fluid_particle_pointers\n");

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

    // Allocate edge index arrays
    edges.edge_pointers_left = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    edges.edge_pointers_right = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    // Allocate out of bound index arrays
    out_of_bounds.oob_pointer_indicies_left = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.oob_pointer_indicies_right = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.vacant_indicies = malloc(2*out_of_bounds.max_oob_particles * sizeof(int));

    printf("bytes allocated: %lu\n", total_bytes);

    // Initialize particles
    initParticles(fluid_particle_pointers, fluid_particles, &water_volume_global, start_x,
		  number_particles_x, &edges, max_fluid_particles_local, spacing_particle, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, params.number_fluid_particles_local, params.tunable_params.smoothing_radius);

    // Send intiial paramaters to render node
    tunable_parameters *null_tunable_param = NULL;
    int *null_recvcnts = NULL;
    int *null_displs = NULL;
    MPI_Gatherv(&params.tunable_params, 1, TunableParamtype, null_tunable_param, null_recvcnts, null_displs, TunableParamtype, 0, MPI_COMM_WORLD);

    // Initialize RGB Light if present
    #if defined BLINK1
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

    int sub_step = 0; // substep range from 0 to < steps_per_frame

    // Main simulation loop
    while(1) {

        // Initialize velocities
        apply_gravity(fluid_particle_pointers, &params);

        // Advance to predicted position and set OOB particles
        predict_positions(fluid_particle_pointers, &boundary_global, &params);

        // Make sure that async send to render node is complete
        if(sub_step == 0)
        {
            if(coords_req != MPI_REQUEST_NULL)
	        MPI_Wait(&coords_req, MPI_STATUS_IGNORE);
        }

        #if defined BLINK1
        char previously_active = params.tunable_params.active;
        #endif

        // Receive updated paramaters from render nodes
        if(sub_step == steps_per_frame-1)
            MPI_Scatterv(null_tunable_param, 0, null_displs, TunableParamtype, &params.tunable_params, 1, TunableParamtype, 0,  MPI_COMM_WORLD);

        #if defined BLINK1
        // If recently added to computation turn light to light state color
        // If recently taken out of computation turn light to white
        char currently_active = params.tunable_params.active;
        if (!previously_active && currently_active)
            rgb_light_reset(&light_state);
        else if (!currently_active && previously_active)
            rgb_light_white(&light_state);
        #endif

        if(params.tunable_params.kill_sim)
            break;

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(fluid_particle_pointers, fluid_particles, &out_of_bounds, &boundary_global, &params);

        // Hash the non halo regions
        hash_fluid(fluid_particle_pointers, &neighbor_grid, &params);

         // Exchange halo particles
        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);
        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

        // Add the halo particles to neighbor buckets
        hash_halo(fluid_particle_pointers, &neighbor_grid, &params);

        int solve_iterations = 4;
        int si;
        for(si=0; si<solve_iterations; si++)
        {
            compute_densities(fluid_particle_pointers, neighbors, &params);

            calculate_lambda(fluid_particle_pointers, neighbors, &params);
            // Generally not needed it appears but included for correctness of parallel algorithm
            update_halo_lambdas(fluid_particle_pointers, &edges, &params);

            update_dp(fluid_particle_pointers, neighbors, &params);

            update_dp_positions(fluid_particle_pointers, &boundary_global, &params);
            // Generally not needed it appears but included for correctness of parallel algorithm
            update_halo_positions(fluid_particle_pointers, &edges, &params);
        }

        // update velocity
        updateVelocities(fluid_particle_pointers, &edges, &boundary_global, &params);

//        vorticity_confinement(fluid_particle_pointers, neighbors, &params);

        XSPH_viscosity(fluid_particle_pointers, neighbors, &params);

        // Update position
        update_positions(fluid_particle_pointers, &params);

        // Pack fluid particle coordinates
        // This sends results as short in pixel coordinates
        if(sub_step == steps_per_frame-1)
        {
            for(i=0; i<params.number_fluid_particles_local; i++) {
                p = fluid_particle_pointers[i];
                fluid_particle_coords[i*2] = (2.0f*p->x/boundary_global.max_x - 1.0f) * SHRT_MAX; // convert to short using full range
                fluid_particle_coords[(i*2)+1] = (2.0f*p->y/boundary_global.max_y - 1.0f) * SHRT_MAX; // convert to short using full range
            }
            // Async send fluid particle coordinates to render node
            MPI_Isend(fluid_particle_coords, 2*params.number_fluid_particles_local, MPI_SHORT, 0, 17, MPI_COMM_WORLD, &coords_req);
        }

        if(sub_step == steps_per_frame-1)
            sub_step = 0;
        else
	    sub_step++;

    }

    #if defined BLINK1
        shutdown_rgb_light(&light_state);
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

// Smoothing kernels

// (h^2 - r^2)^3 normalized in 2D
float W(float r, float h)
{
    if(r > h)
        return 0.0f;

    float C = 4.0f/(M_PI*h*h*h*h*h*h*h*h);
    float W = C*(h*h-r*r)*(h*h-r*r)*(h*h-r*r);
    return W;
}

// Gradient (h-r)^3 normalized in 2D
float del_W(float r, float h)
{
    if(r > h)
        return 0.0f;

    float eps  = 0.000001f;
    float coef = -30.0f/M_PI;
    float C = coef/(h*h*h*h*h * (r+eps));
    float del_W = C*(h-r)*(h-r);
    return del_W;
}

void vorticity_confinement(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;
    float epsilon = 20.01f;
    float dt = params->tunable_params.time_step;

    float x_diff, y_diff, vx_diff, vy_diff, r_mag, dw, dw_x, dw_y, part_vort_z, vort_z, eta_x, eta_y, eta_mag, N_x, N_y;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        vort_z = 0.0f;
        eta_x = 0.0f;
        eta_y = 0.0f;
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = n->fluid_neighbors[j];
            x_diff = p->x_star - q->x_star;
            y_diff = p->y_star - q->y_star;
            vx_diff = q->v_x - p->v_x;
            vy_diff = q->v_y - p->v_y;
            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff);

            dw = del_W(r_mag, params->tunable_params.smoothing_radius);
            dw_x = dw*x_diff;
            dw_y = dw*y_diff;

            part_vort_z = vx_diff*dw_y - vy_diff*dw_x;
            vort_z += part_vort_z;

            if(x_diff<0.0000001 || y_diff<0.0000001)
                continue;

            eta_x += abs(part_vort_z)/x_diff;
            eta_y += abs(part_vort_z)/y_diff;            
        }

        eta_mag = sqrt(eta_x*eta_x + eta_y*eta_y);

        if(eta_mag<0.0000001)
            continue;

        N_x = eta_x / eta_mag;
        N_y = eta_y / eta_mag;

        p->v_x += epsilon*dt*N_y*vort_z;
        p->v_y -= epsilon*dt*N_x*vort_z;
    }
}

void XSPH_viscosity(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;
    float c = 0.1f;

    float x_diff, y_diff, vx_diff, vy_diff, r_mag, w;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        float partial_sum_x = 0.0f;
        float partial_sum_y = 0.0f;
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = n->fluid_neighbors[j];
            x_diff = p->x_star - q->x_star;
            y_diff = p->y_star - q->y_star;
            vx_diff = q->v_x - p->v_x;
            vy_diff = q->v_y - p->v_y;
            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff);
            w = W(r_mag, params->tunable_params.smoothing_radius);
            partial_sum_x += vx_diff * w;
            partial_sum_y += vy_diff * w;
        }
        partial_sum_x *= c;
        partial_sum_y *= c;
        p->v_x += partial_sum_x;
        p->v_y += partial_sum_y;
    }

}

void compute_densities(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        p->density = 0.0f;

        // Own contribution to density
        calculate_density(p,p,params->tunable_params.smoothing_radius, params->particle_mass);
        // Neighbor contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = n->fluid_neighbors[j];
            calculate_density(p, q, params->tunable_params.smoothing_radius, params->particle_mass);
        }
    }

}

// This should go into the hash, perhaps with the viscocity?
void apply_gravity(fluid_particle **fluid_particle_pointers, param *params)
{
    int i;
    fluid_particle *p;
    float dt = params->tunable_params.time_step;
    float g = -params->tunable_params.g;

    for(i=0; i<(params->number_fluid_particles_local); i++) {
        p = fluid_particle_pointers[i];
        p->v_y += g*dt;
     }
}

void update_dp_positions(fluid_particle **fluid_particle_pointers, AABB_t* boundary_global,  param *params)
{
     int i;
     fluid_particle *p;

     for(i=0; i<(params->number_fluid_particles_local); i++) {
        p = fluid_particle_pointers[i];
        p->x_star += p->dp_x;
        p->y_star += p->dp_y;

	// Enforce boundary conditions
        boundaryConditions(p, boundary_global, params);
    }    
}

void update_positions(fluid_particle **fluid_particle_pointers, param *params)
{
     int i;
     fluid_particle *p;

     for(i=0; i<(params->number_fluid_particles_local); i++) {
        p = fluid_particle_pointers[i];
        p->x = p->x_star;
        p->y = p->y_star;
    }    
}

void calculate_lambda(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;

    for(i=0; i<params->number_fluid_particles_local; i++)
    { 
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        float Ci = p->density/params->tunable_params.rest_density - 1.0f;

        float sum_C, x_diff, y_diff, r_mag, grad, grad_x, grad_y;

        sum_C = 0.0f;
        grad_x = 0.0f;
        grad_y = 0.0f;
        // Add k = i contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = n->fluid_neighbors[j];
            x_diff = p->x_star - q->x_star;
            y_diff = p->y_star - q->y_star;
            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff);
            grad = del_W(r_mag, params->tunable_params.smoothing_radius);
            grad_x += grad*x_diff ;
            grad_y += grad*y_diff;
           }
           sum_C += grad_x*grad_x + grad_y*grad_y;

        // Add k =j contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = n->fluid_neighbors[j];
            x_diff = p->x_star - q->x_star;
            y_diff = p->y_star - q->y_star;
            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff);
            grad = del_W(r_mag, params->tunable_params.smoothing_radius);
            grad_x = grad*x_diff ;
            grad_y = grad*y_diff;
            sum_C += (grad_x*grad_x + grad_y*grad_y);
        }

        sum_C *= (1.0f/params->tunable_params.rest_density*params->tunable_params.rest_density);  

        float epsilon = 1.0f;
        p->lambda = -Ci/(sum_C + epsilon);
    }
}

void update_dp(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{

    fluid_particle *p, *q;
    neighbor *n;
    float x_diff, y_diff, dp, r_mag;

    int i,j;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        float dp_x = 0.0f;
        float dp_y = 0.0f;
        float s_corr;
        float k = 0.1f;
        float dq = 0.3f*params->tunable_params.smoothing_radius;
        float Wdq = W(dq, params->tunable_params.smoothing_radius);

        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = n->fluid_neighbors[j];
            x_diff = p->x_star - q->x_star;
            y_diff = p->y_star - q->y_star;
            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff);
            s_corr = -k*(powf(W(r_mag, params->tunable_params.smoothing_radius)/Wdq, 4.0f));
            dp = (p->lambda + q->lambda + s_corr)*del_W(r_mag, params->tunable_params.smoothing_radius);
            dp_x += dp*x_diff;
            dp_y += dp*y_diff;
        }
        p->dp_x = dp_x/params->tunable_params.rest_density;
        p->dp_y = dp_y/params->tunable_params.rest_density;
    }   
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


// Predict position
void predict_positions(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params)
{
    int i;
    fluid_particle *p;
    float dt = params->tunable_params.time_step;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
	p->x_star = p->x  + (p->v_x * dt);
        p->y_star = p->y + (p->v_y * dt);

	// Enforce boundary conditions
        boundaryConditions(p, boundary_global, params);
    }
}

// Calculate the density contribution of p on q and q on p
// r is passed in as this function is called in the hash which must also calculate r
void calculate_density(fluid_particle *p, fluid_particle *q, float h, float mass)
{
    float x_diff, y_diff, r_mag;
    x_diff = p->x_star - q->x_star;
    y_diff = p->y_star - q->y_star;
    r_mag = sqrt(x_diff*x_diff + y_diff*y_diff);
    if(r_mag <= h)
        p->density += mass*W(r_mag, h);
}

void checkVelocity(float *v_x, float *v_y)
{
    const float v_max = 20.0f;

    if(*v_x > v_max)
        *v_x = v_max;
    else if(*v_x < -v_max)
        *v_x = -v_max;
    if(*v_y > v_max)
        *v_y = v_max;
    else if(*v_y < -v_max)
        *v_y = -v_max;
}

void updateVelocity(fluid_particle *p, param *params)
{
    float dt = params->tunable_params.time_step;
    float v_x, v_y;

    v_x = (p->x_star-p->x)/dt;
    v_y = (p->y_star-p->y)/dt;

    checkVelocity(&v_x, &v_y);

    p->v_x = v_x;
    p->v_y = v_y;
}

// Update particle position and check boundary
void updateVelocities(fluid_particle **fluid_particle_pointers, edge_t *edges, AABB_t *boundary_global, param *params)
{
    int i;
    fluid_particle *p;

    // Update local and halo particles, update halo so that XSPH visc. is correct
    for(i=0; i<params->number_fluid_particles_local + params->number_halo_particles; i++) {
        p = fluid_particle_pointers[i];
        updateVelocity(p, params);
    }
}

// Assume AABB with min point being axis origin
void boundaryConditions(fluid_particle *p, AABB_t *boundary, param *params)
{

    float center_x = params->tunable_params.mover_center_x;
    float center_y = params->tunable_params.mover_center_y;

    // Boundary condition for sphere mover
    // Sphere width == height
    float radius = params->tunable_params.mover_width*0.5f;
    float norm_x; 
    float norm_y;

    // Test if inside of circle
    float d;
    float d2 = (p->x_star - center_x)*(p->x_star - center_x) + (p->y_star - center_y)*(p->y_star - center_y);
    if(d2 <= radius*radius && d2 > 0.0f) {
        d = sqrt(d2);
        norm_x = (center_x-p->x_star)/d;
        norm_y = (center_y-p->y_star)/d;
	    
        // With no collision impulse we can handle penetration here
        float pen_dist = radius - d;
        p->x_star -= pen_dist * norm_x;
        p->y_star -= pen_dist * norm_y;
    }

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(p->x_star < boundary->min_x) {
        p->x_star = boundary->min_x;
    }
    else if(p->x_star > boundary->max_x){
        p->x_star = boundary->max_x-0.001f;
    }
    if(p->y_star <  boundary->min_y) {
        p->y_star = boundary->min_y;
    }
    else if(p->y_star > boundary->max_y){
        p->y_star = boundary->max_y-0.001f;
    }
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
        fluid_particle_pointers[i]->v_x = 0.0f;
        fluid_particle_pointers[i]->v_y = 0.0f;
    }
}
