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
#include "communication.h"
#include "fluid.h"
#include <stddef.h>
#include "mpi.h"

// Externally defined in header
MPI_Comm MPI_COMM_COMPUTE;
MPI_Group group_world;
MPI_Group group_compute;
MPI_Group group_render;
// Externally defined in header
MPI_Datatype TunableParamtype;
MPI_Datatype LeftEdgetype;
MPI_Datatype RightEdgetype;

// Rank 0 is the render node, the rest are compute nodes
// This will create appropriate MPI communicators
void create_communicators()
{
    // Extract group handle
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);    

    // Add ranks > 0 to group_compute
    int exclude_rank = 0;
    MPI_Group_excl(group_world, 1, &exclude_rank, &group_compute);

    // Create render group
    int include_rank = 0;
    MPI_Group_incl(group_world, 1, &include_rank, &group_render);

    // Create communicator from group_compute
    MPI_Comm_create(MPI_COMM_WORLD, group_compute, &MPI_COMM_COMPUTE);
}

void create_MPI_types()
{
    MPI_Datatype types[30];
    MPI_Aint disps[30];
    int blocklens[30];
    int i; 

    // Create param type
    for(i=0; i<13; i++) types[i] = MPI_FLOAT;
    types[13] = MPI_CHAR;
    types[14] = MPI_CHAR;
    for (i=0; i<15; i++) blocklens[i] = 1;
    // Get displacement of each struct member
    disps[0] = offsetof( tunable_parameters_t, rest_density );
    disps[1] = offsetof( tunable_parameters_t, smoothing_radius );
    disps[2] = offsetof( tunable_parameters_t, g );
    disps[3] = offsetof( tunable_parameters_t, k );
    disps[4] = offsetof( tunable_parameters_t, dq );
    disps[5] = offsetof( tunable_parameters_t, c );
    disps[6] = offsetof( tunable_parameters_t, time_step );
    disps[7] = offsetof( tunable_parameters_t, proc_start );
    disps[8] = offsetof( tunable_parameters_t, proc_end );
    disps[9] = offsetof( tunable_parameters_t, mover_center_x );
    disps[10] = offsetof( tunable_parameters_t, mover_center_y );
    disps[11] = offsetof( tunable_parameters_t, mover_center_z );
    disps[12] = offsetof( tunable_parameters_t, mover_radius );
    disps[13] = offsetof( tunable_parameters_t, kill_sim );
    disps[14] = offsetof( tunable_parameters_t, active );

    // Commit type
    MPI_Type_create_struct( 15, blocklens, disps, types, &TunableParamtype );
    MPI_Type_commit( &TunableParamtype );
}

void free_MPI_types()
{
    MPI_Type_free(&TunableParamtype);

    MPI_Group_free(&group_world);
    MPI_Group_free(&group_compute);
}

void update_halo_lambdas(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;

    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    int num_moving_left = edges->number_edge_particles_left;
    int num_moving_right = edges->number_edge_particles_right;

    int num_from_left = params->number_halo_particles_left;
    int num_from_right = params->number_halo_particles_right;

    // Allocate send/recv buffers
    // Could combine left/right into single malloc...
    float *send_lambdas_left = (float*)malloc(sizeof(float)*num_moving_left);
    float *send_lambdas_right = (float*)malloc(sizeof(float)*num_moving_right);

    float *recv_lambdas_left = (float*)malloc(sizeof(float)*num_from_left);
    float *recv_lambdas_right = (float*)malloc(sizeof(float)*num_from_right);

    // Pack local halo lambdas
    for(i=0; i<num_moving_left; i++) {
        p_index = edges->edge_indices_left[i];
        send_lambdas_left[i] = fluid_particles->lambda[p_index];
    }
    for(i=0; i<num_moving_right; i++) {
        p_index = edges->edge_indices_right[i];
        send_lambdas_right[i] = fluid_particles->lambda[p_index];
    }

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Could do async to perhaps increase performance
    // Send lambdas to right and receive from left
    int tag = 784;
    MPI_Sendrecv(send_lambdas_right, num_moving_right, MPI_FLOAT, proc_to_right, tag,
                 recv_lambdas_left, num_from_left, MPI_FLOAT, proc_to_left, tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    // Send lambdas to left and receive from right
    tag = 456;
    MPI_Sendrecv(send_lambdas_left, num_moving_left, MPI_FLOAT, proc_to_left, tag,
                 recv_lambdas_right, num_from_right, MPI_FLOAT, proc_to_right,tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    // Unpack halo particle lambdas
    for(i=0; i<num_from_left; i++) {
        p_index = fluid_particle_indices[params->number_fluid_particles_local + i];;
        fluid_particles->lambda[p_index] = recv_lambdas_left[i];
    }
    for(i=0; i<num_from_right; i++) {
        p_index = fluid_particle_indices[params->number_fluid_particles_local + num_from_left + i];
        fluid_particles->lambda[p_index] = recv_lambdas_right[i];
    }

    // Cleanup memory
    free(send_lambdas_left);
    free(send_lambdas_right);
    free(recv_lambdas_left);
    free(recv_lambdas_right);
}

void update_halo_positions(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;

    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    int num_moving_left = 3*edges->number_edge_particles_left;
    int num_moving_right = 3*edges->number_edge_particles_right;

    int num_from_left = 3*params->number_halo_particles_left;
    int num_from_right = 3*params->number_halo_particles_right;

    // Allocate send/recv buffers
    // Could combine left/right into single malloc...
    float *send_positions_left = (float*)malloc(sizeof(float)*num_moving_left);
    float *send_positions_right = (float*)malloc(sizeof(float)*num_moving_right);

    float *recv_positions_left = (float*)malloc(sizeof(float)*num_from_left);
    float *recv_positions_right = (float*)malloc(sizeof(float)*num_from_right);

    // Pack local edge positions
    for(i=0; i<num_moving_left; i+=3) {
        p_index = edges->edge_indices_left[i/3];
        send_positions_left[i] = fluid_particles->x_star[p_index];
        send_positions_left[i+1] = fluid_particles->y_star[p_index];
        send_positions_left[i+2] = fluid_particles->z_star[p_index];
    }
    for(i=0; i<num_moving_right; i+=3) {
        p_index = edges->edge_indices_right[i/3];
        send_positions_right[i] = fluid_particles->x_star[p_index];
        send_positions_right[i+1] = fluid_particles->y_star[p_index];
        send_positions_right[i+2] = fluid_particles->z_star[p_index];
    }

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Could do async to perhaps increase performance
    // Send positions to right and receive from left
    int tag = 874;
    MPI_Sendrecv(send_positions_right, num_moving_right, MPI_FLOAT, proc_to_right, tag,
                 recv_positions_left, num_from_left, MPI_FLOAT, proc_to_left, tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    // Send positions to left and receive from right
    tag = 546;
    MPI_Sendrecv(send_positions_left, num_moving_left, MPI_FLOAT, proc_to_left, tag,
                 recv_positions_right, num_from_right, MPI_FLOAT, proc_to_right,tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    // Unpack halo particle positions
    for(i=0; i<num_from_left; i+=3) {
        p_index = fluid_particle_indices[params->number_fluid_particles_local + i/3];;
        fluid_particles->x_star[p_index] = recv_positions_left[i];
        fluid_particles->y_star[p_index] = recv_positions_left[i+1];
        fluid_particles->z_star[p_index] = recv_positions_left[i+2];
    }
    for(i=0; i<num_from_right; i+=3) {
        p_index = fluid_particle_indices[params->number_fluid_particles_local + num_from_left/3 + i/3];
        fluid_particles->x_star[p_index] = recv_positions_right[i];
        fluid_particles->y_star[p_index] = recv_positions_right[i+1];
        fluid_particles->z_star[p_index] = recv_positions_right[i+2];
    }

    // Cleanup memory
    free(send_positions_left);
    free(send_positions_right);
    free(recv_positions_left);
    free(recv_positions_right);
}

// Pack particle struct float components into contiguous memory
void pack_halo_components(float *left_send, float *right_send, fluid_sim_t *fluid_sim)
{
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;

    // Append halos, fluid particle struct has 14 float components
    int i, p_index;
    for (i=0; i<edges->number_edge_particles_left; i++) {
        p_index = edges->edge_indices_left[i];
        left_send[i*14]     = fluid_particles->x_star[p_index];
        left_send[i*14 + 1] = fluid_particles->y_star[p_index];
        left_send[i*14 + 2] = fluid_particles->z_star[p_index];
        left_send[i*14 + 3] = fluid_particles->x[p_index];
        left_send[i*14 + 4] = fluid_particles->y[p_index];
        left_send[i*14 + 5] = fluid_particles->z[p_index];
        left_send[i*14 + 6] = fluid_particles->v_x[p_index];
        left_send[i*14 + 7] = fluid_particles->v_y[p_index];
        left_send[i*14 + 8] = fluid_particles->v_z[p_index];
        left_send[i*14 + 9] = fluid_particles->dp_x[p_index];
        left_send[i*14 + 10] = fluid_particles->dp_y[p_index];
        left_send[i*14 + 11] = fluid_particles->dp_z[p_index];
        left_send[i*14 + 12] = fluid_particles->density[p_index];
        left_send[i*14 + 13] = fluid_particles->lambda[p_index];
    }
    for (i=0; i<edges->number_edge_particles_right; i++) {
        p_index = edges->edge_indices_right[i];
        right_send[i*14]     = fluid_particles->x_star[p_index];
        right_send[i*14 + 1] = fluid_particles->y_star[p_index];
        right_send[i*14 + 2] = fluid_particles->z_star[p_index];
        right_send[i*14 + 3] = fluid_particles->x[p_index];
        right_send[i*14 + 4] = fluid_particles->y[p_index];
        right_send[i*14 + 5] = fluid_particles->z[p_index];
        right_send[i*14 + 6] = fluid_particles->v_x[p_index];
        right_send[i*14 + 7] = fluid_particles->v_y[p_index];
        right_send[i*14 + 8] = fluid_particles->v_z[p_index];
        right_send[i*14 + 9] = fluid_particles->dp_x[p_index];
        right_send[i*14 + 10] = fluid_particles->dp_y[p_index];
        right_send[i*14 + 11] = fluid_particles->dp_z[p_index];
        right_send[i*14 + 12] = fluid_particles->density[p_index];
        right_send[i*14 + 13] = fluid_particles->lambda[p_index];
    }
}

// Unpack halo components
void unpack_halo_components(float *packed_recv_left, float *packed_recv_right, fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    param_t *params = fluid_sim->params;
    edge_t *edges = fluid_sim->edges;

    int i;
    uint p_index;

    // Unpack halo particles from left rank first
    for(i=0; i<params->number_halo_particles_left; i++)
    {
        p_index = params->max_fluid_particle_index + 1 + i; // "Global" index
        fluid_particles->x_star[p_index]  = packed_recv_left[i*14];
        fluid_particles->y_star[p_index]  = packed_recv_left[i*14 + 1];
        fluid_particles->z_star[p_index]  = packed_recv_left[i*14 + 2];
        fluid_particles->x[p_index]       = packed_recv_left[i*14 + 3];
        fluid_particles->y[p_index]       = packed_recv_left[i*14 + 4];
        fluid_particles->z[p_index]       = packed_recv_left[i*14 + 5];
        fluid_particles->v_x[p_index]     = packed_recv_left[i*14 + 6];
        fluid_particles->v_y[p_index]     = packed_recv_left[i*14 + 7];
        fluid_particles->v_z[p_index]     = packed_recv_left[i*14 + 8];
        fluid_particles->dp_x[p_index]    = packed_recv_left[i*14 + 9];
        fluid_particles->dp_y[p_index]    = packed_recv_left[i*14 + 10];
        fluid_particles->dp_z[p_index]    = packed_recv_left[i*14 + 11];
        fluid_particles->density[p_index] = packed_recv_left[i*14 + 12];
        fluid_particles->lambda[p_index]  = packed_recv_left[i*14 + 13];
        fluid_particles->id[p_index]      = params->number_fluid_particles_local + i;
        fluid_particle_indices[params->number_fluid_particles_local+i] = p_index;
    }

    // Unpack halo particles from right rank second
    for(i=0; i<params->number_halo_particles_right; i++)
    {
        p_index = params->max_fluid_particle_index + 1 + params->number_halo_particles_left + i; // "Global" index
        fluid_particles->x_star[p_index]  = packed_recv_right[i*14];
        fluid_particles->y_star[p_index]  = packed_recv_right[i*14 + 1];
        fluid_particles->z_star[p_index]  = packed_recv_right[i*14 + 2];
        fluid_particles->x[p_index]       = packed_recv_right[i*14 + 3];
        fluid_particles->y[p_index]       = packed_recv_right[i*14 + 4];
        fluid_particles->z[p_index]       = packed_recv_right[i*14 + 5];
        fluid_particles->v_x[p_index]     = packed_recv_right[i*14 + 6];
        fluid_particles->v_y[p_index]     = packed_recv_right[i*14 + 7];
        fluid_particles->v_z[p_index]     = packed_recv_right[i*14 + 8];
        fluid_particles->dp_x[p_index]    = packed_recv_right[i*14 + 9];
        fluid_particles->dp_y[p_index]    = packed_recv_right[i*14 + 10];
        fluid_particles->dp_z[p_index]    = packed_recv_right[i*14 + 11];
        fluid_particles->density[p_index] = packed_recv_right[i*14 + 12];
        fluid_particles->lambda[p_index]  = packed_recv_right[i*14 + 13];
        fluid_particles->id[p_index]      = params->number_fluid_particles_local + params->number_halo_particles_left + i;
        fluid_particle_indices[params->number_fluid_particles_local + params->number_halo_particles_left + i] = p_index;
   }
}

void halo_exchange(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;

    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;
    float h = params->tunable_params.smoothing_radius;

    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    // Set edge particle indicies and update number
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p_index = fluid_particle_indices[i];
        if (fluid_particles->x[p_index] - params->tunable_params.proc_start <= h)
            edges->edge_indices_left[edges->number_edge_particles_left++] = p_index;
        else if (params->tunable_params.proc_end - fluid_particles->x[p_index] <= h)
            edges->edge_indices_right[edges->number_edge_particles_right++] = p_index;
    }

    int num_moving_left = edges->number_edge_particles_left;
    int num_moving_right = edges->number_edge_particles_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    debug_print("rank %d, halo: will send %d to left, %d to right\n", rank, num_moving_left, num_moving_right);

    // Get number of halo particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;

    // Send number of particles to right and receive from left
    int tag = 3217;
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);
    // Send number of particles to left and receive from right
    tag = 8425;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    debug_print("rank %d, halo: will recv %d from left, %d from right\n", rank, num_from_left, num_from_right);

    // Allocate memory for packed arrays
    int num_components = 14;
    float *packed_send_left = (float*)malloc(num_components * num_moving_left * sizeof(float));
    float *packed_recv_left = (float*)malloc(num_components * num_from_left * sizeof(float));
    float *packed_send_right = (float*)malloc(num_components * num_moving_right * sizeof(float));
    float *packed_recv_right = (float*)malloc(num_components * num_from_right * sizeof(float));

    // Pack halo particle struct components to send
    pack_halo_components(packed_send_left, packed_send_right, fluid_sim);

    debug_print("rank %d, prams->max_fluid_particle_index: %d\n", rank,  params->max_fluid_particle_index);
    debug_print("rank %d, halo: send %d to left, %d to right\n", rank, num_moving_left, num_moving_right);

    tag = 4312;
    // Send packed particles to right and receive from left
    MPI_Sendrecv(packed_send_right, num_components*num_moving_right, MPI_FLOAT, proc_to_right, tag, 
                 packed_recv_left,  num_components*num_from_left   , MPI_FLOAT, proc_to_left, tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);
    tag = 5177;
    // Send packed particles to left and receive from right
    MPI_Sendrecv(packed_send_left,  num_components*num_moving_left, MPI_FLOAT, proc_to_left, tag,
                 packed_recv_right, num_components*num_from_right,  MPI_FLOAT, proc_to_right,tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    // Need to automatically add rank to debug print
    debug_print("rank %d halo: recv %d from left, %d from right\n",rank, num_from_left,num_from_right);

    // Update params struct with halo values
    int total_received = num_from_left + num_from_right;
    params->number_halo_particles = total_received;
    params->number_halo_particles_left = num_from_left;
    params->number_halo_particles_right = num_from_right;

    // Unpack halo components from left and right
    unpack_halo_components(packed_recv_left, packed_recv_right, fluid_sim);

    // Free memory
    free(packed_send_left);
    free(packed_recv_left);
    free(packed_send_right);
    free(packed_recv_right);
}

// Pack out of bounds particle components
void pack_oob_components(float *left_send, float *right_send, fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    oob_t *oob = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

    // Append halos, fluid particle struct has 14 float components
    int i, p_index;
    for (i=0; i<oob->number_oob_particles_left; i++) {
        p_index = fluid_particle_indices[oob->oob_index_indices_left[i]];
        left_send[i*14]     = fluid_particles->x_star[p_index];
        left_send[i*14 + 1] = fluid_particles->y_star[p_index];
        left_send[i*14 + 2] = fluid_particles->z_star[p_index];
        left_send[i*14 + 3] = fluid_particles->x[p_index];
        left_send[i*14 + 4] = fluid_particles->y[p_index];
        left_send[i*14 + 5] = fluid_particles->z[p_index];
        left_send[i*14 + 6] = fluid_particles->v_x[p_index];
        left_send[i*14 + 7] = fluid_particles->v_y[p_index];
        left_send[i*14 + 8] = fluid_particles->v_z[p_index];
        left_send[i*14 + 9] = fluid_particles->dp_x[p_index];
        left_send[i*14 + 10] = fluid_particles->dp_y[p_index];
        left_send[i*14 + 11] = fluid_particles->dp_z[p_index];
        left_send[i*14 + 12] = fluid_particles->density[p_index];
        left_send[i*14 + 13] = fluid_particles->lambda[p_index];
        // Invalidate index entry as particle is now gone
        fluid_particle_indices[oob->oob_index_indices_left[i]] = ((uint)-1);

        // Add index to array of vacancies
        oob->vacant_indices[oob->number_vacancies++] = p_index;

        params->number_fluid_particles_local--;
    }
    for (i=0; i<oob->number_oob_particles_right; i++) {
        p_index = fluid_particle_indices[oob->oob_index_indices_right[i]];
        right_send[i*14]     = fluid_particles->x_star[p_index];
        right_send[i*14 + 1] = fluid_particles->y_star[p_index];
        right_send[i*14 + 2] = fluid_particles->z_star[p_index];
        right_send[i*14 + 3] = fluid_particles->x[p_index];
        right_send[i*14 + 4] = fluid_particles->y[p_index];
        right_send[i*14 + 5] = fluid_particles->z[p_index];
        right_send[i*14 + 6] = fluid_particles->v_x[p_index];
        right_send[i*14 + 7] = fluid_particles->v_y[p_index];
        right_send[i*14 + 8] = fluid_particles->v_z[p_index];
        right_send[i*14 + 9] = fluid_particles->dp_x[p_index];
        right_send[i*14 + 10] = fluid_particles->dp_y[p_index];
        right_send[i*14 + 11] = fluid_particles->dp_z[p_index];
        right_send[i*14 + 12] = fluid_particles->density[p_index];
        right_send[i*14 + 13] = fluid_particles->lambda[p_index];
        // Invalidate index entry as particle is now gone
        fluid_particle_indices[oob->oob_index_indices_right[i]] = ((uint)-1);

        // Add index to array of vacancies
        oob->vacant_indices[oob->number_vacancies++] = p_index;

        // Decrement number of fluid particles
        params->number_fluid_particles_local--;
    }
}

// Unpack out of bounds components
void unpack_oob_components(float *packed_recv, int num_recv, fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    oob_t *oob = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

    int i;

    // We must fill in invalid entries in the particle_index array
    // Invalid entries are caused by particles that have just left
    int num_invalid_left = oob->number_oob_particles_left;
    int num_invalid_right = oob->number_oob_particles_right;
    int num_invalid_total = num_invalid_left + num_invalid_right;

    uint p_index;

    // Number of invalidated index entries that have been replaced
    int indices_replaced = 0;

    // Unpack oob particles into vacancies if possible
    for(i=0; i<num_recv; i++)
    {
        // Unpack into vacancies first
        if(oob->number_vacancies > 0) {
            p_index = oob->vacant_indices[oob->number_vacancies-1];
            oob->number_vacancies--;
        }
        else // If no vacancies add onto end of global array
            p_index = ++params->max_fluid_particle_index;

        // Update pointer array
        if(i < num_invalid_total) { // If there are invalid entries in pointer array update those first
            if(indices_replaced < num_invalid_left) 
                fluid_particle_indices[oob->oob_index_indices_left[indices_replaced]] = p_index; 
            else
                fluid_particle_indices[oob->oob_index_indices_right[i-num_invalid_left]] = p_index;

            indices_replaced ++;
        }
        else { // If no invalid entries add to end of index array and increase number of local particles
            fluid_particle_indices[params->number_fluid_particles_local] = p_index;
        }

        // Incriment number of local fluid particles
        params->number_fluid_particles_local++;

        fluid_particles->x_star[p_index]  = packed_recv[i*14];
        fluid_particles->y_star[p_index]  = packed_recv[i*14 + 1];
        fluid_particles->z_star[p_index]  = packed_recv[i*14 + 2];
        fluid_particles->x[p_index]       = packed_recv[i*14 + 3];
        fluid_particles->y[p_index]       = packed_recv[i*14 + 4];
        fluid_particles->z[p_index]       = packed_recv[i*14 + 5];
        fluid_particles->v_x[p_index]     = packed_recv[i*14 + 6];
        fluid_particles->v_y[p_index]     = packed_recv[i*14 + 7];
        fluid_particles->v_z[p_index]     = packed_recv[i*14 + 8];
        fluid_particles->dp_x[p_index]    = packed_recv[i*14 + 9];
        fluid_particles->dp_y[p_index]    = packed_recv[i*14 + 10];
        fluid_particles->dp_z[p_index]    = packed_recv[i*14 + 11];
        fluid_particles->density[p_index] = packed_recv[i*14 + 12];
        fluid_particles->lambda[p_index]  = packed_recv[i*14 + 13];
    }
}

// Transfer particles that are out of node bounds
void transfer_OOB_particles(fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    oob_t *out_of_bounds = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

    int i;
    uint p_index;

    int rank;
    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    int nprocs;
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    int num_moving_left = out_of_bounds->number_oob_particles_left;
    int num_moving_right = out_of_bounds->number_oob_particles_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;
    int tag = 7006;
    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_COMPUTE,MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8278;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_COMPUTE,MPI_STATUS_IGNORE);

    // Allocate memory for packed arrays
    int num_components = 14;
    float *packed_send_left = (float*)malloc(num_components * num_moving_left * sizeof(float));
    float *packed_send_right = (float*)malloc(num_components * num_moving_right * sizeof(float));
    float *packed_recv = (float*)malloc(num_components * (num_from_right+num_from_left) * sizeof(float));

    // Pack OOB particle struct components to send
    pack_oob_components(packed_send_left, packed_send_right, fluid_sim);

    // Send packed particles to right and receive from left
    tag = 2522;
    MPI_Sendrecv(packed_send_right, num_components*num_moving_right, MPI_FLOAT, proc_to_right, tag,
                 packed_recv,  num_components*num_from_left   , MPI_FLOAT, proc_to_left, tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);
    tag = 1165;
    // Send packed particles to left and receive from right
    MPI_Sendrecv(packed_send_left,  num_components*num_moving_left, MPI_FLOAT, proc_to_left, tag,
                 (packed_recv + num_components*num_from_left), num_components*num_from_right,  MPI_FLOAT, proc_to_right,tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    int total_received = num_from_right + num_from_left;
    int total_sent     = num_moving_right + num_moving_left;

    // Unpack components and update vacancies for particles that were just received
    unpack_oob_components(packed_recv, total_received, fluid_sim);

    debug_print("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n", rank, num_moving_left, num_moving_right, num_from_left, num_from_right);
    debug_print("rank %d OOB: num vacant %d\n", rank, out_of_bounds->number_vacancies);

    // Go through all possible fluid particles and remove null entries
    int num_particles = 0;
    int diff = 0;
    // If more particles were sent than recvd we need to check more indices than number_fluid_particles_local
    if(total_sent>total_received)
        diff = total_sent-total_received;
    int check_length = params->number_fluid_particles_local + diff; //The max length in the index array we need to check
    for (i=0; i<check_length; i++) {
        p_index = fluid_particle_indices[i];
        if (p_index != (uint)-1) {
            fluid_particle_indices[num_particles] = p_index;
            fluid_particles->id[p_index] = num_particles;
            num_particles++;
        }
    }

    // Need to add rank to debug_print
    debug_print("rank %d num local: %d\n", rank, num_particles);
    debug_print("rank %d params->num local %d\n", rank, params->number_fluid_particles_local);

    // Free memory
    free(packed_send_left);
    free(packed_send_right);
    free(packed_recv);
}

