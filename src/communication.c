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
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
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
        p = edges->edge_pointers_left[i];
        send_lambdas_left[i] = p->lambda;
    }
    for(i=0; i<num_moving_right; i++) {
        p = edges->edge_pointers_right[i];
        send_lambdas_right[i] = p->lambda;
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
        p = fluid_particle_pointers[params->number_fluid_particles_local + i];;
        p->lambda = recv_lambdas_left[i];
    }
    for(i=0; i<num_from_right; i++) {
        p = fluid_particle_pointers[ params->number_fluid_particles_local + num_from_left + i];
        p->lambda = recv_lambdas_right[i];
    }

    // Cleanup memory
    free(send_lambdas_left);
    free(send_lambdas_right);
    free(recv_lambdas_left);
    free(recv_lambdas_right);
}

void update_halo_positions(fluid_sim_t *fluid_sim)
{
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    int i;
    fluid_particle_t *p;

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
        p = edges->edge_pointers_left[i/2];
        send_positions_left[i] = p->x_star;
        send_positions_left[i+1] = p->y_star;
    }
    for(i=0; i<num_moving_right; i+=2) {
        p = edges->edge_pointers_right[i/2];
        send_positions_right[i] = p->x_star;
        send_positions_right[i+1] = p->y_star;
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
        p = fluid_particle_pointers[params->number_fluid_particles_local + i/2];;
        p->x_star = recv_positions_left[i];
        p->y_star = recv_positions_left[i+1];
    }
    for(i=0; i<num_from_right; i+=2) {
        p = fluid_particle_pointers[ params->number_fluid_particles_local + num_from_left/2 + i/2];
        p->x_star = recv_positions_right[i];
        p->y_star = recv_positions_right[i+1];
    }

    // Cleanup memory
    free(send_positions_left);
    free(send_positions_right);
    free(recv_positions_left);
    free(recv_positions_right);
}

void start_halo_exchange(fluid_sim_t *fluid_sim)
{
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
    fluid_particle_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    int i;
    fluid_particle_t *p;
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
        p = fluid_particle_pointers[i];
        if (p->x - params->tunable_params.node_start_x <= h)
            edges->edge_pointers_left[edges->number_edge_particles_left++] = p;
        else if (params->tunable_params.node_end_x - p->x <= h)
            edges->edge_pointers_right[edges->number_edge_particles_right++] = p;
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

    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8425;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_COMPUTE, MPI_STATUS_IGNORE);

    debug_print("rank %d, halo: will recv %d from left, %d from right\n", rank, num_from_left, num_from_right);

    int *blocklens_left = (int*)malloc(num_moving_left * sizeof(int));
    int *blocklens_right = (int*)malloc(num_moving_right * sizeof(int));
    int *indicies_left = (int*)malloc(num_moving_left * sizeof(int));
    int *indicies_right = (int*)malloc(num_moving_right * sizeof(int));

    // Convert the edge pointer into a particle array index using pointer arithmetic
    int index;
    int edge_global_index; // Should perhaps use ptrdiff_t
    for (i=0; i<num_moving_left; i++) {
        blocklens_left[i] = 1;
        edge_global_index = (int) (edges->edge_pointers_left[i] - fluid_particles);
        indicies_left[i]  = edge_global_index;
    }
    for (i=0; i<num_moving_right; i++) {
        blocklens_right[i] = 1;
        edge_global_index = (int) (edges->edge_pointers_right[i] - fluid_particles);
        indicies_right[i]  = edge_global_index;
    }

    MPI_Type_indexed(num_moving_left,blocklens_left,indicies_left,Particletype,&LeftEdgetype);
    MPI_Type_commit(&LeftEdgetype);
    MPI_Type_indexed(num_moving_right,blocklens_right,indicies_right,Particletype,&RightEdgetype);
    MPI_Type_commit(&RightEdgetype);

    debug_print("rank %d, prams->max_fluid_particle_index: %d\n", rank,  params->max_fluid_particle_index);

    //Index to start receiving halo particles
    int indexToReceiveLeft = params->max_fluid_particle_index + 1;
    int indexToReceiveRight = indexToReceiveLeft + num_from_left;

    debug_print("rank %d, halo: send %d to left, %d to right, indexToRecieveLeft %d, indexToReceiveRight %d \n", rank, num_moving_left, num_moving_right, indexToReceiveLeft, indexToReceiveRight);

    int tagl = 4312;
    int tagr = 5177;
    // Receive halo from left rank
    MPI_Irecv(&fluid_particles[indexToReceiveLeft], num_from_left, Particletype, proc_to_left,tagl, MPI_COMM_COMPUTE, &edges->reqs[0]);
    // Receive halo from right rank
    MPI_Irecv(&fluid_particles[indexToReceiveRight], num_from_right, Particletype, proc_to_right,tagr, MPI_COMM_COMPUTE, &edges->reqs[1]);
    // Send halo to right rank
    MPI_Isend(fluid_particles,1,RightEdgetype,proc_to_right,tagl,MPI_COMM_COMPUTE, &edges->reqs[2]);
    MPI_Isend(fluid_particles,1,LeftEdgetype,proc_to_left,tagr,MPI_COMM_COMPUTE, &edges->reqs[3]);

    // Free allocated arays
    free(blocklens_left);
    free(blocklens_right);
    free(indicies_left);
    free(indicies_right);
}

void finish_halo_exchange(fluid_sim_t *fluid_sim)
{

    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
    fluid_particle_t *fluid_particles = fluid_sim->fluid_particles;
    edge_t *edges = fluid_sim->edges;
    param_t *params = fluid_sim->params;

    int i;
    // Wait for transfer to complete
    MPI_Status statuses[4];
    MPI_Waitall(4, edges->reqs, statuses);

    int num_received_right = 0;
    int num_received_left = 0;
    MPI_Get_count(&statuses[0], Particletype, &num_received_left);
    MPI_Get_count(&statuses[1], Particletype, &num_received_right);

    int total_received = num_received_left + num_received_right;
    params->number_halo_particles = total_received;
    params->number_halo_particles_left = num_received_left;
    params->number_halo_particles_right = num_received_right;

    // Need to automatically add rank to debug print
    debug_print("halo: recv %d from left, %d from right\n",num_received_left,num_received_right);

    // Update pointer array with new values
    int local_index;
    int global_index;
    for (i=0; i< total_received; i++) {
        local_index = params->number_fluid_particles_local + i;
        global_index = params->max_fluid_particle_index + 1 + i;
        fluid_particle_pointers[local_index] = &fluid_particles[global_index];
        fluid_particle_pointers[local_index]->id = local_index;
    }

    // Free indexed types
    MPI_Type_free(&LeftEdgetype);
    MPI_Type_free(&RightEdgetype);
}

// Transfer particles that are out of node bounds
void transfer_OOB_particles(fluid_sim_t *fluid_sim)
{
    fluid_particle_t **fluid_particle_pointers = fluid_sim->fluid_particle_pointers;
    fluid_particle_t *fluid_particles = fluid_sim->fluid_particles;
    oob_t *out_of_bounds = fluid_sim->out_of_bounds;
    param_t *params = fluid_sim->params;

    int i;
    fluid_particle_t *p;

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

    // Create indexed type to send
    MPI_Datatype LeftMovetype;
    MPI_Datatype RightMovetype;
    int *blocklens_left = malloc(num_moving_left*sizeof(int));
    int *blocklens_right = malloc(num_moving_right*sizeof(int));
    int *indicies_left = malloc(num_moving_left*sizeof(int));
    int *indicies_right = malloc(num_moving_right*sizeof(int));

    // Convert the OOB pointer into a particle array index using pointer arithmetic
    int index;
    int oob_global_index;
    for (i=0; i<num_moving_left; i++) {
        blocklens_left[i] = 1;
        index = out_of_bounds->oob_pointer_indicies_left[i];
        oob_global_index = (int) (fluid_particle_pointers[index] - fluid_particles);
        indicies_left[i]  = oob_global_index;
    }
    for (i=0; i<num_moving_right; i++) {
        blocklens_right[i] = 1;
        index = out_of_bounds->oob_pointer_indicies_right[i];
        oob_global_index = (int) (fluid_particle_pointers[index] - fluid_particles);
        indicies_right[i]  = oob_global_index;
    }

    MPI_Type_indexed(num_moving_left,blocklens_left,indicies_left,Particletype,&LeftMovetype);
    MPI_Type_commit(&LeftMovetype);
    MPI_Type_indexed(num_moving_right,blocklens_right,indicies_right,Particletype,&RightMovetype);
    MPI_Type_commit(&RightMovetype);

    // Create indexed type to recv
    MPI_Datatype LeftRecvtype;
    MPI_Datatype RightRecvtype;

    int total_recv = num_from_left + num_from_right;
    int *blocklens_recv = malloc(total_recv*sizeof(int));
    int *indicies_recv = malloc(total_recv*sizeof(int));
    // Go through vacancies, starting at end, to create receive indicies
    // Go through reverse so it's easy to set update vacancies below
    for (i=0; i<out_of_bounds->number_vacancies && i<total_recv; i++) {
        blocklens_recv[i] = 1;
        index = out_of_bounds->number_vacancies-1-i;
        indicies_recv[i] = out_of_bounds->vacant_indicies[index];
    }
    int replaced = i;
    int remaining = total_recv - replaced;
    for (i=0; i<remaining; i++) {
        blocklens_recv[replaced+i] = 1;
        indicies_recv[replaced+i] = params->max_fluid_particle_index + 1 + i;
    }

    MPI_Type_indexed(num_from_left, blocklens_recv, indicies_recv, Particletype, &LeftRecvtype);
    MPI_Type_commit(&LeftRecvtype);
    MPI_Type_indexed(num_from_right, &blocklens_recv[num_from_left],&indicies_recv[num_from_left], Particletype, &RightRecvtype);
    MPI_Type_commit(&RightRecvtype); 

    MPI_Status status;
    MPI_Request request;

    // Send oob particles to right processor receive oob particles from right processor
    int num_received_left = 0;
    int num_received_right = 0;

    // Sending to right, recv from left
    tag = 2522;
    MPI_Sendrecv(fluid_particles,1,RightMovetype,proc_to_right,tag,fluid_particles,1,LeftRecvtype,proc_to_left,tag,MPI_COMM_COMPUTE,&status);
    MPI_Get_count(&status, Particletype, &num_received_left);
    // Sending to left, recv from right
    tag = 1165;
    MPI_Sendrecv(fluid_particles,1,LeftMovetype,proc_to_left,tag,fluid_particles,1,RightRecvtype,proc_to_right,tag,MPI_COMM_COMPUTE,&status);
    MPI_Get_count(&status, Particletype, &num_received_right);

    int total_sent = num_moving_left + num_moving_right;
    int total_received = num_received_right + num_received_left;

    debug_print("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n", rank, num_moving_left, num_moving_right, num_received_left, num_received_right);

    // Update maximum particle index if neccessary
    int max_received_index = total_received?total_received-1:0;// If non received don't access indicies_recv[-1]...
    if (total_received && indicies_recv[max_received_index] > params->max_fluid_particle_index) {
        debug_print("rank %d increasing max index from: %d", rank, params->max_fluid_particle_index);
        params->max_fluid_particle_index = indicies_recv[max_received_index];
        debug_print("to %d\n", params->max_fluid_particle_index);
    }

    // Update vacancy total for particles received
    if (total_received < out_of_bounds->number_vacancies )
        out_of_bounds->number_vacancies -= total_received;
    else
        out_of_bounds->number_vacancies = 0;

    // Update vacancy indicies and total for particles sent
    // Set sent particle pointer to received particle location or to NULL
    int oob_pointer_index;
    int recv_replaced = 0; // received particles that have replaced leaving particles
    for (i=0; i<num_moving_left; i++) {
        // Index of particle pointer that has left
        oob_pointer_index = out_of_bounds->oob_pointer_indicies_left[i];
        // Index of particle that has left
        oob_global_index = (int) (fluid_particle_pointers[oob_pointer_index] - fluid_particles);
        out_of_bounds->vacant_indicies[out_of_bounds->number_vacancies] = oob_global_index;
        // Incriment the number of vacancies
        out_of_bounds->number_vacancies++;

        // Set pointer from removed particle to a recvd particle ,if any, else NULL
        if(recv_replaced < total_received) {
            fluid_particle_pointers[oob_pointer_index] = &fluid_particles[indicies_recv[recv_replaced]]; 
            recv_replaced++;
        }
        else
            fluid_particle_pointers[oob_pointer_index] = NULL;
    }
    for (i=0; i<num_moving_right; i++) {
        // Index of particle pointer that has left
        oob_pointer_index = out_of_bounds->oob_pointer_indicies_right[i];
        // Index of particle that has left
        oob_global_index = (int) (fluid_particle_pointers[oob_pointer_index] - fluid_particles);
        // Set new vacancy index
        out_of_bounds->vacant_indicies[out_of_bounds->number_vacancies] = oob_global_index;
        // Incriment the number of vacancies
        out_of_bounds->number_vacancies++;

        // Set pointer from removed particle to a recvd particle ,if any, else NULL
        if(recv_replaced < total_received) {
            fluid_particle_pointers[oob_pointer_index] = &fluid_particles[indicies_recv[recv_replaced]];
            recv_replaced++;
        }
        else
            fluid_particle_pointers[oob_pointer_index] = NULL;
    }

    debug_print("rank %d OOB: num vacant %d\n", rank, out_of_bounds->number_vacancies);

    // If more particles are received than sent add to end of pointer array
    remaining = total_received - recv_replaced;
    int pointer_index;
    int global_index;
    int max_fluid_pointers = params->number_fluid_particles_local;
    for(i=0; i<remaining; i++)
    {
        pointer_index = params->number_fluid_particles_local + i;
        global_index = indicies_recv[recv_replaced + i];
        fluid_particle_pointers[pointer_index] = &fluid_particles[global_index];
        max_fluid_pointers++;
    }

    // Update particle pointer array
    // Go through all possible fluid particles and remove null entries
    int num_particles = 0;
    for (i=0; i<max_fluid_pointers; i++) {
        p = fluid_particle_pointers[i];
        if (p != NULL) {
            fluid_particle_pointers[num_particles] = p;
            fluid_particle_pointers[num_particles]->id = num_particles;
            num_particles++;
        }
    }

    params->number_fluid_particles_local = num_particles;

    // Need to add rank to debug_print
    debug_print("num local: %d\n", num_particles);

    // Free indexed types
    MPI_Type_free(&LeftRecvtype);
    MPI_Type_free(&RightRecvtype);
    MPI_Type_free(&LeftMovetype);
    MPI_Type_free(&RightMovetype);

    free(blocklens_left);
    free(blocklens_right);
    free(indicies_left);
    free(indicies_right);
    free(blocklens_recv);
    free(indicies_recv);
}

