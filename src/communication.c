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
#include "mpi.h"
#include "communication.h"
#include "fluid.h"
#include <stddef.h>

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

/*
    // Create fluid particle type;
    for (i=0; i<10; i++) types[i] = MPI_FLOAT;
    types[10] = MPI_INT;
    for (i=0; i<11; i++) blocklens[i] = 1;
    // Get displacement of each struct member
    disps[0] = offsetof( fluid_particle_t, x_star);
    disps[1] = offsetof( fluid_particle_t, y_star);
    disps[2] = offsetof( fluid_particle_t, x);
    disps[3] = offsetof( fluid_particle_t, y);
    disps[4] = offsetof( fluid_particle_t, v_x);
    disps[5] = offsetof( fluid_particle_t, v_y);
    disps[6] = offsetof( fluid_particle_t, density);
    disps[7]   = offsetof( fluid_particle_t, dp_x);
    disps[8]   = offsetof( fluid_particle_t, dp_y);
    disps[9]   = offsetof( fluid_particle_t, lambda);
    disps[10] = offsetof( fluid_particle_t, id);
    // Commit type
    MPI_Type_create_struct( 11, blocklens, disps, types, &Particletype );
    MPI_Type_commit( &Particletype );
*/
    // Create param type
    for(i=0; i<15; i++) types[i] = MPI_FLOAT;
    types[15] = MPI_CHAR;
    types[16] = MPI_CHAR;
    for (i=0; i<17; i++) blocklens[i] = 1;
    // Get displacement of each struct member
    disps[0] = offsetof( tunable_parameters_t, rest_density );
    disps[1] = offsetof( tunable_parameters_t, smoothing_radius );
    disps[2] = offsetof( tunable_parameters_t, g );
    disps[3] = offsetof( tunable_parameters_t, k );
    disps[4] = offsetof( tunable_parameters_t, k_near );
    disps[5] = offsetof( tunable_parameters_t, k_spring );
    disps[6] = offsetof( tunable_parameters_t, sigma );
    disps[7] = offsetof( tunable_parameters_t, beta );
    disps[8] = offsetof( tunable_parameters_t, time_step );
    disps[9] = offsetof( tunable_parameters_t, node_start_x );
    disps[10] = offsetof( tunable_parameters_t, node_end_x );
    disps[11] = offsetof( tunable_parameters_t, mover_center_x );
    disps[12] = offsetof( tunable_parameters_t, mover_center_y );
    disps[13] = offsetof( tunable_parameters_t, mover_width );
    disps[14] = offsetof( tunable_parameters_t, mover_height );
    disps[15] = offsetof( tunable_parameters_t, kill_sim );
    disps[16] = offsetof( tunable_parameters_t, active );

    // Commit type
    MPI_Type_create_struct( 17, blocklens, disps, types, &TunableParamtype );
    MPI_Type_commit( &TunableParamtype );
}

void free_MPI_types()
{
    MPI_Type_free(&Particletype);
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
    fluid_particle_t *p;

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
    float *send_lambdas_left = malloc(sizeof(float)*num_moving_left);
    float *send_lambdas_right = malloc(sizeof(float)*num_moving_right);

    float *recv_lambdas_left = malloc(sizeof(float)*num_from_left);
    float *recv_lambdas_right = malloc(sizeof(float)*num_from_right);

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

    int num_moving_left = 2*edges->number_edge_particles_left;
    int num_moving_right = 2*edges->number_edge_particles_right;

    int num_from_left = 2*params->number_halo_particles_left;
    int num_from_right = 2*params->number_halo_particles_right;

    // Allocate send/recv buffers
    // Could combine left/right into single malloc...
    float *send_positions_left = malloc(sizeof(float)*num_moving_left);
    float *send_positions_right = malloc(sizeof(float)*num_moving_right);

    float *recv_positions_left = malloc(sizeof(float)*num_from_left);
    float *recv_positions_right = malloc(sizeof(float)*num_from_right);

    // Pack local edge positions
    for(i=0; i<num_moving_left; i+=2) {
        p_index = edges->edge_indices_left[i/2];
        send_positions_left[i] = fluid_particles->x_star[p_index];
        send_positions_left[i+1] = fluid_particles->y_star[p_index];
    }
    for(i=0; i<num_moving_right; i+=2) {
        p_index = edges->edge_pointers_right[i/2];
        send_positions_right[i] = fluid_particles->x_star[p_index];
        send_positions_right[i+1] = fluid_particles->y_star[p_index];
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
    for(i=0; i<num_from_left; i+=2) {
        p_index = fluid_particle_indices[params->number_fluid_particles_local + i/2];;
        fluid_particles->x_star[p_index] = recv_positions_left[i];
        fluid_particles->y_star[p_index] = recv_positions_left[i+1];
    }
    for(i=0; i<num_from_right; i+=2) {
        p = fluid_particle_indices[params->number_fluid_particles_local + num_from_left/2 + i/2];
        fluid_particles->x_star[p_index] = recv_positions_right[i];
        fluid_particles->y_star[p_index] = recv_positions_right[i+1];
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
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;

    // Append halos, fluid particle struct has 10 float components
    int i, p_index;
    for (i=0; i<num_moving_left; i++) {
        p_index = edges->edge_indices_left[i];
        left_send[i*10]     = fluid_particles->x_star[p_index];
        left_send[i*10 + 1] = fluid_particles->y_star[p_index];
        left_send[i*10 + 2] = fluid_particles->x[p_index];
        left_send[i*10 + 3] = fluid_particles->y[p_index];
        left_send[i*10 + 4] = fluid_particles->v_x[p_index];
        left_send[i*10 + 5] = fluid_particles->v_y[p_index];
        left_send[i*10 + 6] = fluid_particles->dp_x[p_index];
        left_send[i*10 + 7] = fluid_particles->dp_y[p_index];
        left_send[i*10 + 8] = fluid_particles->density[p_index];
        left_send[i*10 + 9] = fluid_particles->lambda[p_index];
    }
    for (i=0; i<num_moving_right; i++) {
        p_index = edges->edge_indices_right[i];
        right_send[i*10]     = fluid_particles->x_star[p_index];
        right_send[i*10 + 1] = fluid_particles->y_star[p_index];
        right_send[i*10 + 2] = fluid_particles->x[p_index];
        right_send[i*10 + 3] = fluid_particles->y[p_index];
        right_send[i*10 + 4] = fluid_particles->v_x[p_index];
        right_send[i*10 + 5] = fluid_particles->v_y[p_index];
        right_send[i*10 + 6] = fluid_particles->dp_x[p_index];
        right_send[i*10 + 7] = fluid_particles->dp_y[p_index];
        right_send[i*10 + 8] = fluid_particles->density[p_index];
        right_send[i*10 + 9] = fluid_particles->lambda[p_index];
    }
}

// Unpack halo components
void unpack_halo_components(packed_recv_left, packed_recv_right, edges, fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;

    int i;
    uint p_index;
    // Unpack halo particles from left rank first
    for(i=0; i<params->number_halo_particles_left; i++)
    {
        p_index = params->max_fluid_particle_index + 1 + i; // "Global" index
        fluid_particles->x_star[p_index]  = packed_recv_left[i*10];
        fluid_particles->y_star[p_index]  = packed_recv_left[i*10 + 1];
        fluid_particles->x[p_index]       = packed_recv_left[i*10 + 2];
        fluid_particles->y[p_index]       = packed_recv_left[i*10 + 3];
        fluid_particles->v_x[p_index]     = packed_recv_left[i*10 + 4];
        fluid_particles->v_y[p_index]     = packed_recv_left[i*10 + 5];
        fluid_particles->dp_x[p_index]    = packed_recv_left[i*10 + 6];
        fluid_particles->dp_y[p_index]    = packed_recv_left[i*10 + 7];
        fluid_particles->density[p_index] = packed_recv_left[i*10 + 8];
        fluid_particles->lambda[p_index]  = packed_recv_left[i*10 + 9];
        fluid_particles->id[p_index]      = params->number_fluid_particles_local + i;
        fluid_particle_indices[params->number_fluid_particles_local+i] = p_index;
    }

    // Unpack halo particles from right rank second
    for(i=0; i<params->number_halo_particles_right; i++)
    {
        p_index = params->max_fluid_particle_index + 1 + params->number_halo_particles_left + i; // "Global" index
        fluid_particles->x_star[p_index]  = packed_recv_right[i*10];
        fluid_particles->y_star[p_index]  = packed_recv_right[i*10 + 1];
        fluid_particles->x[p_index]       = packed_recv_right[i*10 + 2];
        fluid_particles->y[p_index]       = packed_recv_right[i*10 + 3];
        fluid_particles->v_x[p_index]     = packed_recv_right[i*10 + 4];
        fluid_particles->v_y[p_index]     = packed_recv_right[i*10 + 5];
        fluid_particles->dp_x[p_index]    = packed_recv_right[i*10 + 6];
        fluid_particles->dp_y[p_index]    = packed_recv_right[i*10 + 7];
        fluid_particles->density[p_index] = packed_recv_right[i*10 + 8];
        fluid_particles->lambda[p_index]  = packed_recv_right[i*10 + 9];
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
        if (fluid_particles->x[p_index] - params->tunable_params.node_start_x <= h)
            edges->edge_indices_left[edges->number_edge_particles_left++] = p_index;
        else if (params->tunable_params.node_end_x - fluid_particles->x[p_index] <= h)
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
    int tag = 3217;

    // Send number of particles to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);
    // Send number of particles to left and receive from right
    tag = 8425;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    debug_print("rank %d, halo: will recv %d from left, %d from right\n", rank, num_from_left, num_from_right);

    // Allocate memory for packed arrays
    int num_components = 10;
    float *packed_send_left = (float*)malloc(num_components * num_moving_left * sizeof(float));
    float *packed_recv_left = (float*)malloc(num_components * num_from_left * sizeof(float));
    float *packed_send_right = (float*)malloc(num_components * num_moving_right * sizeof(float));
    float *packed_recv_right = (float*)malloc(num_components * num_from_right * sizeof(float));

    // Pack halo particle struct components to send
    pack_halo_components(packed_left_send, packed_right_send, fluid_sim);

    debug_print("rank %d, prams->max_fluid_particle_index: %d\n", rank,  params->max_fluid_particle_index);
    debug_print("rank %d, halo: send %d to left, %d to right, indexToRecieveLeft %d, indexToReceiveRight %d \n", rank, num_moving_left, num_moving_right, indexToReceiveLeft, indexToReceiveRight);

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
    debug_print("halo: recv %d from left, %d from right\n",num_received_left,num_received_right);

    // Update params struct with halo values
    int total_received = num_from_left + num_from_right;
    params->number_halo_particles = total_received;
    params->number_halo_particles_left = num_from_left;
    params->number_halo_particles_right = num_from_right;

    // Unpack halo components from left and right
    unpack_halo_components(packed_left_recv, packed_right_recv, fluid_sim);

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
    oob_t *oob = fluid_sim->oob;

    // Append halos, fluid particle struct has 10 float components
    int i, p_index;
    for (i=0; i<params->number_oob_particles_left; i++) {
        p_index = oob->oob_indices_left[i];
        left_send[i*10]     = fluid_particles->x_star[p_index];
        left_send[i*10 + 1] = fluid_particles->y_star[p_index];
        left_send[i*10 + 2] = fluid_particles->x[p_index];
        left_send[i*10 + 3] = fluid_particles->y[p_index];
        left_send[i*10 + 4] = fluid_particles->v_x[p_index];
        left_send[i*10 + 5] = fluid_particles->v_y[p_index];
        left_send[i*10 + 6] = fluid_particles->dp_x[p_index];
        left_send[i*10 + 7] = fluid_particles->dp_y[p_index];
        left_send[i*10 + 8] = fluid_particles->density[p_index];
        left_send[i*10 + 9] = fluid_particles->lambda[p_index];
    }
    for (i=0; i<params->number_oob_particles_right; i++) {
        p_index = oob->oob_indices_right[i];
        right_send[i*10]     = fluid_particles->x_star[p_index];
        right_send[i*10 + 1] = fluid_particles->y_star[p_index];
        right_send[i*10 + 2] = fluid_particles->x[p_index];
        right_send[i*10 + 3] = fluid_particles->y[p_index];
        right_send[i*10 + 4] = fluid_particles->v_x[p_index];
        right_send[i*10 + 5] = fluid_particles->v_y[p_index];
        right_send[i*10 + 6] = fluid_particles->dp_x[p_index];
        right_send[i*10 + 7] = fluid_particles->dp_y[p_index];
        right_send[i*10 + 8] = fluid_particles->density[p_index];
        right_send[i*10 + 9] = fluid_particles->lambda[p_index];
    }
}

// Unpack out of bounds components
void unpack_oob_components(float *packed_recv, int num_recv, edge_t *edges, fluid_sim_t *fluid_sim)
{
    uint *fluid_particle_indices = fluid_sim->fluid_particle_indices;
    fluid_particles_t *fluid_particles = fluid_sim->fluid_particles;
    oob_t *oob = fluid_sim->oob;

    int i;
    uint p_index;
    // Unpack oob particles
    for(i=0; i<num_recv; i++)
    {
        // Unpack into vacancies first
        if(oob->number_vacancies > 0)
            p_index = fluid_particle_indices[oob->vancant_indices[oob->number_vacancies--]];
        else
            p_index = ++params->max_fluid_particle_index;

        fluid_particles->x_star[p_index]  = packed_recv[i*10];
        fluid_particles->y_star[p_index]  = packed_recv[i*10 + 1];
        fluid_particles->x[p_index]       = packed_recv[i*10 + 2];
        fluid_particles->y[p_index]       = packed_recv[i*10 + 3];
        fluid_particles->v_x[p_index]     = packed_recv[i*10 + 4];
        fluid_particles->v_y[p_index]     = packed_recv[i*10 + 5];
        fluid_particles->dp_x[p_index]    = packed_recv[i*10 + 6];
        fluid_particles->dp_y[p_index]    = packed_recv[i*10 + 7];
        fluid_particles->density[p_index] = packed_recv[i*10 + 8];
        fluid_particles->lambda[p_index]  = packed_recv[i*10 + 9];
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
    int num_components = 10;
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
                 (packed_recv + num_from_left), num_components*num_from_right,  MPI_FLOAT, proc_to_right,tag,
                 MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    int total_sent = num_moving_left + num_moving_right;
    int total_received = num_from_right + num_from_left;

    // Update vacancies to include particles that were just sent
    for (i=0; i<num_moving_left; i++) {
        // Index in particle index array of particle that has left
        oob_index = (int) out_of_bounds->oob_index_indices_left[i];
        // Add index to array of vacancies
        out_of_bounds->vacant_indices[out_of_bounds->number_vacancies] = oob_index;
        // Incriment the number of vacancies
        out_of_bounds->number_vacancies++;
        // Invalidate index as particle is now gone
        fluid_particle_indices[oob_index] = ((uint)-1);
    }
    for (i=0; i<num_moving_right; i++) {
        // Index in particle index array of particle that has left
        oob_index = (int) out_of_bounds->oob_index_indices_right[i];
        // Add index to array of vacancies
        out_of_bounds->vacant_indices[out_of_bounds->number_vacancies] = oob_index;
        // Incriment the number of vacancies
        out_of_bounds->number_vacancies++;
        // Invalidate index as particle is now gone
        fluid_particle_indices[oob_index] = ((uint)-1);
    }

    // Unpack components and update vacancies for particles that were just received
    unpack_oob_components(packed_recv, total_received, fluid_sim);

    debug_print("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n", rank, num_moving_left, num_moving_right, num_from_left, num_from_right);
    debug_print("rank %d OOB: num vacant %d\n", rank, out_of_bounds->number_vacancies);

    // If more particles received than sent update shit

    // Update particle pointer array
    // Go through all possible fluid particles and remove null entries
    int num_particles = 0;
    for (i=0; i<max_fluid_pointers; i++) {
        p = fluid_particle_pointers[i];
        if (p != (uint)-1) {
            fluid_particle_pointers[num_particles] = p;
            fluid_particle_pointers[num_particles]->id = num_particles;
            num_particles++;
        }
    }

    params->number_fluid_particles_local = num_particles;

    // Need to add rank to debug_print
    debug_print("num local: %d\n", num_particles);
}

