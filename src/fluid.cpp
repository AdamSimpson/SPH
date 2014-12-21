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

#include "hash_sort.h"
#include "renderer.h"
#include "setup.h"
#include "fluid.h"
#include "communication.h"

#ifdef BLINK1
#include "blink1_light.h"
#endif

void start_simulation()
{
    unsigned int i;

    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    printf("compute rank: %d, num compute procs: %d \n",rank, nprocs);

    // Struct to hold all simulation values
    fluid_sim_t fluid_sim;

    // Allocate structs used in simulation
    alloc_sim_structs(&fluid_sim);

    // Initialize simulation parameters
    // Additionally set fluid and  boundaries
    // Requires MPI_Bcast to send screen aspect ratio
    init_params(&fluid_sim);

    // Partition problem, allocate memory, and initialize particles
    // Requires MPI_Send for particle counts
    alloc_and_init_sim(&fluid_sim);

    // Send initial parameters to render node and initialize light
    // Requires MPI_Gatherv
    sync_initial_params(&fluid_sim);

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

    // Setup dummy values for MPI
    float *null_float = NULL;

    MPI_Request coords_req = MPI_REQUEST_NULL;

    tunable_parameters_t *null_tunable_param = NULL;
    int *null_recvcnts = NULL;
    int *null_displs = NULL;
    uint p_index;

    int sub_step = 0; // substep range from 0 to < steps_per_frame

    // Main simulation loop
    while(1) {

        // Initialize velocities
        apply_gravity(&fluid_sim);

        // Advance to predicted position and set OOB particles
        predict_positions(&fluid_sim);

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
        if(sub_step == fluid_sim.params->steps_per_frame-1)
            MPI_Scatterv(null_tunable_param, 0, null_displs, TunableParamtype, &fluid_sim.params->tunable_params, 1, TunableParamtype, 0,  MPI_COMM_WORLD);

        #if defined BLINK1
        // If recently added to computation turn light to light state color
        // If recently taken out of computation turn light to white
        char currently_active = fluid_sim.params->tunable_params.active;
        if (!previously_active && currently_active)
            rgb_light_reset(&light_state);
        else if (!currently_active && previously_active)
            rgb_light_white(&light_state);
        #endif

        if(fluid_sim.params->tunable_params.kill_sim)
            break;

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(&fluid_sim);

         // Exchange halo particles
        halo_exchange(&fluid_sim);

        // Hash particles, sort, fill particle neighbors
        find_all_neighbors(&fluid_sim);

        int solve_iterations = 4;
        int si;
        for(si=0; si<solve_iterations; si++)
        {
            compute_densities(&fluid_sim);

            calculate_lambda(&fluid_sim);
            // Generally not needed it appears but included for correctness of parallel algorithm
            update_halo_lambdas(&fluid_sim);

            update_dp(&fluid_sim);

            update_dp_positions(&fluid_sim);
            // Generally not needed it appears but included for correctness of parallel algorithm
            update_halo_positions(&fluid_sim);
        }

        // update velocity
        update_velocities(&fluid_sim);

//        vorticity_confinement(fluid_particle_pointers, neighbors, &params);

        XSPH_viscosity(&fluid_sim);

        // Update position
        update_positions(&fluid_sim);

        // Pack fluid particle coordinates
        // This sends results as short to speed up transfer a bit
        // We need a range of -1 to 1 so we assume max_x is largest value
        // and so we normalize x to -1 to 1 and the rest are scaled by the same value
        if(sub_step == fluid_sim.params->steps_per_frame-1)
        {
            for(i=0; i<fluid_sim.params->number_fluid_particles_local; i++) {
                p_index = fluid_sim.fluid_particle_indices[i];
                fluid_sim.fluid_particle_coords[i*3] = (2.0f*fluid_sim.fluid_particles->x[p_index]
                                                       /fluid_sim.boundary_global->max_x - 1.0f) * SHRT_MAX; // convert to short using full range
                fluid_sim.fluid_particle_coords[(i*3)+1] = (2.0f*fluid_sim.fluid_particles->y[p_index]
                                                           /fluid_sim.boundary_global->max_x - 1.0f) * SHRT_MAX; // convert to short using full range
                fluid_sim.fluid_particle_coords[(i*3)+2] = (2.0f*fluid_sim.fluid_particles->z[p_index]
                                                           /fluid_sim.boundary_global->max_x - 1.0f) * SHRT_MAX; // convert to short using full range
            }
            // Async send fluid particle coordinates to render node
            MPI_Isend(fluid_sim.fluid_particle_coords, 3*fluid_sim.params->number_fluid_particles_local, MPI_SHORT, 0, 17, MPI_COMM_WORLD, &coords_req);
        }

        if(sub_step == fluid_sim.params->steps_per_frame-1)
            sub_step = 0;
        else
	    sub_step++;

    }

    #if defined BLINK1
        shutdown_rgb_light(&light_state);
    #endif

    // Free main sim memory
    free_sim_memory(&fluid_sim);

    // Cleanup structs
    free_sim_structs(&fluid_sim);
}

// Smoothing kernels

// (h^2 - r^2)^3 normalized in 3D (poly6)
float W(float r, float h)
{
    if(r > h)
        return 0.0f;

    double C = 315.0/(64.0*M_PI*pow((double)h,9.0));
    float W = C*pow((double)(h*h-r*r), 3.0);
    return W;
}

// Gradient (h-r)^3 normalized in 3D (Spikey)
float del_W(float r, float h)
{
    float C = -45.0/(M_PI * pow((double)h, 6.0));
    float del_W = C*(h-r)*(h-r);
    return del_W;
}

/*
void vorticity_confinement(fluid_sim_t *fluid_sim)
{
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers; 
    neighbor_t *neighbors = fluid_sim->neighbor_grid->neighbors;
    param_t *params = fluid_sim->params;

    int i,j;
    fluid_particle_t *p, *q;
    neighbor_t *n;
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
*/

void XSPH_viscosity(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    neighbor_t *neighbors = fluid_sim->neighbor_grid->neighbors;
    param_t *params = fluid_sim->params;

    int i,j;
    uint p_index, q_index;
    neighbor_t *n;
    float c = params->tunable_params.c;

    float x_diff, y_diff, z_diff, vx_diff, vy_diff, vz_diff, r_mag, w;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p_index = fluid_particle_indices[i];
        n = &neighbors[i];

        float partial_sum_x = 0.0f;
        float partial_sum_y = 0.0f;
        float partial_sum_z = 0.0f;
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles->x_star[p_index] - fluid_particles->x_star[q_index];
            y_diff = fluid_particles->y_star[p_index] - fluid_particles->y_star[q_index];
            z_diff = fluid_particles->z_star[p_index] - fluid_particles->z_star[q_index];

            vx_diff = fluid_particles->v_x[q_index] - fluid_particles->v_x[p_index];
            vy_diff = fluid_particles->v_y[q_index] - fluid_particles->v_y[p_index];
            vz_diff = fluid_particles->v_z[q_index] - fluid_particles->v_z[p_index];

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff * z_diff*z_diff);
            w = W(r_mag, params->tunable_params.smoothing_radius);
            partial_sum_x += vx_diff * w;
            partial_sum_y += vy_diff * w;
            partial_sum_z += vz_diff * w;
        }
        partial_sum_x *= c;
        partial_sum_y *= c;
        partial_sum_z *= c;

        fluid_particles->v_x[p_index] += partial_sum_x;
        fluid_particles->v_y[p_index] += partial_sum_y;
        fluid_particles->v_z[p_index] += partial_sum_z;
    }
}

void compute_densities(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    neighbor_t *neighbors = fluid_sim->neighbor_grid->neighbors;
    param_t *params = fluid_sim->params;

    int i,j;
    uint p_index, q_index;
    neighbor_t *n;
    float h = params->tunable_params.smoothing_radius;
    float mass = 1.0f;// Should just remove mass

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p_index = fluid_particle_indices[i];
        n = &neighbors[i];

        float x_diff, y_diff, z_diff, r_mag, density;
        density = 0.0f;

        // Own contribution to density
        density += mass*W(0.0f, h);

        // Neighbor contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles->x_star[p_index] - fluid_particles->x_star[q_index];
            y_diff = fluid_particles->y_star[p_index] - fluid_particles->y_star[q_index];
            z_diff = fluid_particles->z_star[p_index] - fluid_particles->z_star[q_index];

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            if(r_mag <= h)
                density += mass*W(r_mag, h);
        }

        // Update particle density
        fluid_particles->density[p_index] = density;
    }
}

void apply_gravity(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;
    float dt = params->tunable_params.time_step;
    float g = -params->tunable_params.g;

    for(i=0; i<(params->number_fluid_particles_local); i++) {
        p_index = fluid_particle_indices[i];
        fluid_particles->v_y[p_index] += g*dt;
     }
}

void update_dp_positions(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;

    for(i=0; i<(params->number_fluid_particles_local); i++) {
        p_index = fluid_particle_indices[i];
        fluid_particles->x_star[p_index] += fluid_particles->dp_x[p_index];
        fluid_particles->y_star[p_index] += fluid_particles->dp_y[p_index];
        fluid_particles->z_star[p_index] += fluid_particles->dp_z[p_index];

	// Enforce boundary conditions
        boundary_conditions(p_index, fluid_sim);
    }    
}

void update_positions(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    param_t *params = fluid_sim->params;

     int i;
     uint p_index;

     for(i=0; i<(params->number_fluid_particles_local); i++) {
        p_index = fluid_particle_indices[i];
        fluid_particles->x[p_index] = fluid_particles->x_star[p_index];
        fluid_particles->y[p_index] = fluid_particles->y_star[p_index];
        fluid_particles->z[p_index] = fluid_particles->z_star[p_index];
    }    
}

void calculate_lambda(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    neighbor_t *neighbors = fluid_sim->neighbor_grid->neighbors;
    param_t *params = fluid_sim->params;

    int i,j;
    uint p_index, q_index;
    neighbor_t *n;

    for(i=0; i<params->number_fluid_particles_local; i++)
    { 
        p_index = fluid_particle_indices[i];
        n = &neighbors[i];

        float Ci = fluid_particles->density[p_index]/params->tunable_params.rest_density - 1.0f;

        float sum_C, x_diff, y_diff, z_diff, r_mag, grad, grad_x, grad_y, grad_z;

        sum_C = 0.0f;
        grad_x = 0.0f;
        grad_y = 0.0f;
        grad_z = 0.0f;
        // Add k = i contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles->x_star[p_index] - fluid_particles->x_star[q_index];
            y_diff = fluid_particles->y_star[p_index] - fluid_particles->y_star[q_index];
            z_diff = fluid_particles->z_star[p_index] - fluid_particles->z_star[q_index];

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            grad = del_W(r_mag, params->tunable_params.smoothing_radius);
            grad_x += grad*x_diff;
            grad_y += grad*y_diff;
            grad_z += grad*z_diff;
           }
           sum_C += grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;

        // Add k =j contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles->x_star[p_index] - fluid_particles->x_star[q_index];
            y_diff = fluid_particles->y_star[p_index] - fluid_particles->y_star[q_index];
            z_diff = fluid_particles->z_star[p_index] - fluid_particles->z_star[q_index];

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            grad = del_W(r_mag, params->tunable_params.smoothing_radius);
            grad_x = grad*x_diff ;
            grad_y = grad*y_diff;
            grad_z = grad*z_diff;
            sum_C += (grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
        }

        sum_C *= (1.0f/(params->tunable_params.rest_density*params->tunable_params.rest_density));  

        float epsilon = 1.0f;
        fluid_particles->lambda[p_index] = -Ci/(sum_C + epsilon);
    }
}

void update_dp(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    neighbor_t *neighbors = fluid_sim->neighbor_grid->neighbors;
    param_t *params = fluid_sim->params;

    uint p_index, q_index;
    neighbor_t *n;
    float x_diff, y_diff, z_diff, dp, r_mag;

    int i,j;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p_index = fluid_particle_indices[i];
        n = &neighbors[i];

        float dp_x = 0.0f;
        float dp_y = 0.0f;
        float dp_z = 0.0f;
        float s_corr;
        float k = params->tunable_params.k;
        float dq = params->tunable_params.dq;
        float Wdq = W(dq, params->tunable_params.smoothing_radius);

        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles->x_star[p_index] - fluid_particles->x_star[q_index];
            y_diff = fluid_particles->y_star[p_index] - fluid_particles->y_star[q_index];
            z_diff = fluid_particles->z_star[p_index] - fluid_particles->z_star[q_index];

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            s_corr = -k*(powf(W(r_mag, params->tunable_params.smoothing_radius)/Wdq, 4.0f));
            dp = (fluid_particles->lambda[p_index] + fluid_particles->lambda[q_index] + s_corr)*del_W(r_mag, params->tunable_params.smoothing_radius);
            dp_x += dp*x_diff;
            dp_y += dp*y_diff;
            dp_z += dp*z_diff;
        }
        fluid_particles->dp_x[p_index] = dp_x/params->tunable_params.rest_density;
        fluid_particles->dp_y[p_index] = dp_y/params->tunable_params.rest_density;
        fluid_particles->dp_z[p_index] = dp_z/params->tunable_params.rest_density;
    }   
}

// Identify out of bounds particles and send them to appropriate rank
void identify_oob_particles(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    oob_t *out_of_bounds = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p_index = fluid_particle_indices[i];

        // Set OOB particle indices and update number
        if (fluid_particles->x[p_index] < params->tunable_params.node_start_x)
            out_of_bounds->oob_index_indices_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (fluid_particles->x[p_index] > params->tunable_params.node_end_x)
            out_of_bounds->oob_index_indices_right[out_of_bounds->number_oob_particles_right++] = i;
    }
 
   // Transfer particles that have left the processor bounds
   transfer_OOB_particles(fluid_sim);
}


// Predict position
void predict_positions(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;
    float dt = params->tunable_params.time_step;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p_index = fluid_particle_indices[i];
	fluid_particles->x_star[p_index] = fluid_particles->x[p_index] + (fluid_particles->v_x[p_index] * dt);
        fluid_particles->y_star[p_index] = fluid_particles->y[p_index] + (fluid_particles->v_y[p_index] * dt);
        fluid_particles->z_star[p_index] = fluid_particles->z[p_index] + (fluid_particles->v_z[p_index] * dt);

	// Enforce boundary conditions
        boundary_conditions(p_index, fluid_sim);
    }
}

void check_velocity(float *v_x, float *v_y, float *v_z)
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
    if(*v_z > v_max)
        *v_z = v_max;
    else if(*v_z < -v_max)
        *v_z = -v_max;

}

// Update particle position and check boundary
void update_velocities(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;

    float dt = params->tunable_params.time_step;
    float v_x, v_y, v_z;

    // Update local and halo particles, update halo so that XSPH visc. is correct
    for(i=0; i<params->number_fluid_particles_local + params->number_halo_particles; i++) {
        p_index = fluid_particle_indices[i];

        v_x = (fluid_particles->x_star[p_index] - fluid_particles->x[p_index])/dt;
        v_y = (fluid_particles->y_star[p_index] - fluid_particles->y[p_index])/dt;
        v_z = (fluid_particles->z_star[p_index] - fluid_particles->z[p_index])/dt;

        check_velocity(&v_x, &v_y, &v_z);

        fluid_particles->v_x[p_index] = v_x;
        fluid_particles->v_y[p_index] = v_y;
        fluid_particles->v_z[p_index] = v_z;
    }
}

// Assume AABB with min point being axis origin
void boundary_conditions(uint p_index, fluid_sim_t *fluid_sim)
{
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    AABB_t *boundary = fluid_sim->boundary_global;
    param_t *params = fluid_sim->params;

    float center_x = params->tunable_params.mover_center_x;
    float center_y = params->tunable_params.mover_center_y;
    float center_z = params->tunable_params.mover_center_z;

    // Boundary condition for sphere mover
    // Sphere width == height
    float radius = params->tunable_params.mover_width*0.5f;
    float norm_x; 
    float norm_y;
    float norm_z;

    // Test if inside of sphere
    float d;
    float d2 = (fluid_particles->x_star[p_index] - center_x)*(fluid_particles->x_star[p_index] - center_x)
             + (fluid_particles->y_star[p_index] - center_y)*(fluid_particles->y_star[p_index] - center_y)
             + (fluid_particles->z_star[p_index] - center_z)*(fluid_particles->z_star[p_index] - center_z);

    if(d2 <= radius*radius && d2 > 0.0f) {
        d = sqrt(d2);
        norm_x = (center_x - fluid_particles->x_star[p_index])/d;
        norm_y = (center_y - fluid_particles->y_star[p_index])/d;
        norm_z = (center_z - fluid_particles->z_star[p_index])/d;
	    
        // With no collision impulse we can handle penetration here
        float pen_dist = radius - d;
        fluid_particles->x_star[p_index] -= pen_dist * norm_x;
        fluid_particles->y_star[p_index] -= pen_dist * norm_y;
        fluid_particles->z_star[p_index] -= pen_dist * norm_z;
    }

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(fluid_particles->x_star[p_index]  < boundary->min_x) {
        fluid_particles->x_star[p_index]  = boundary->min_x;
    }
    else if(fluid_particles->x_star[p_index]  > boundary->max_x){
        fluid_particles->x_star[p_index]  = boundary->max_x-0.001f;
    }
    if(fluid_particles->y_star[p_index]  <  boundary->min_y) {
        fluid_particles->y_star[p_index]  = boundary->min_y;
    }
    else if(fluid_particles->y_star[p_index]  > boundary->max_y){
        fluid_particles->y_star[p_index]  = boundary->max_y-0.001f;
    }
    if(fluid_particles->z_star[p_index]  <  boundary->min_z) {
        fluid_particles->z_star[p_index]  = boundary->min_z;
    }
    else if(fluid_particles->z_star[p_index]  > boundary->max_z){
        fluid_particles->z_star[p_index]  = boundary->max_z-0.001f;
    }
}
