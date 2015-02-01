#include <stdio.h>
#include "mpi.h"
#include "communication.h"
#include "fluid.h"
#include <stddef.h>
#include <string.h>

void createMpiTypes()
{
    //Create fluid particle type
    MPI_Datatype types[16];
    int i;
    for (i=0; i<15; i++) types[i] = MPI_DOUBLE;
    types[15] = MPI_UB;
    int blocklens[16];
    for (i=0; i<16; i++) blocklens[i] = 1;
    MPI_Aint disps[16];
    // Get displacement of each struct member
    disps[0] = offsetof( fluid_particle, x);
    disps[1] = offsetof( fluid_particle, y);
    disps[2] = offsetof( fluid_particle, z);
    disps[3] = offsetof( fluid_particle, v_x);
    disps[4] = offsetof( fluid_particle, v_y);
    disps[5] = offsetof( fluid_particle, v_z);
    disps[6] = offsetof( fluid_particle, v_half_x);
    disps[7] = offsetof( fluid_particle, v_half_y);
    disps[8] = offsetof( fluid_particle, v_half_z);
    disps[9] = offsetof( fluid_particle, a_x);
    disps[10] = offsetof( fluid_particle, a_y);
    disps[11] = offsetof( fluid_particle, a_z);
    disps[12] = offsetof( fluid_particle, density);
    disps[13] = offsetof( fluid_particle, pressure);
    disps[14] = offsetof( fluid_particle, id);
    disps[15] = sizeof(fluid_particle);
    // Commit type
    MPI_Type_create_struct( 16, blocklens, disps, types, &Particletype );
    MPI_Type_commit( &Particletype );
}

void freeMpiTypes()
{
    MPI_Type_free(&Particletype);
}

void startHaloExchange(fluid_particle *fluid_particles,  edge *edges, param *params)
{
    int i;
    int rank = params->rank;
    int nprocs = params->nprocs;

    fluid_particle *p;
    double h = params->smoothing_radius;

    // Set edge particle indices and update number
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = &fluid_particles[i];
        if (p->x - params->node_start_x <= h)
            edges->edge_indices_left[edges->number_edge_particles_left++] = i;
        else if (params->node_end_x - p->x <= h)
            edges->edge_indices_right[edges->number_edge_particles_right++] = i;
    }

    int num_moving_left = edges->number_edge_particles_left;
    int num_moving_right = edges->number_edge_particles_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of halo particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;
    int tag = 3217;
    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8425;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Allocate send buffers
    int total_send = num_moving_left + num_moving_right;
    fluid_particle *send_buffer = malloc(total_send*sizeof(fluid_particle));
    fluid_particle *sendl_buffer = send_buffer;
    fluid_particle *sendr_buffer = send_buffer + num_moving_left;

    // Pack send buffers
    for(i=0; i<edges->number_edge_particles_left; i++)
        sendl_buffer[i] = fluid_particles[edges->edge_indices_left[i]];
    for(i=0; i<edges->number_edge_particles_right; i++)
        sendr_buffer[i] = fluid_particles[edges->edge_indices_right[i]];

    //Index to start receiving halo particles
    int indexToReceiveLeft = params->number_fluid_particles_local;
    int indexToReceiveRight = indexToReceiveLeft + num_from_left;

    int tagl = 4312;
    int tagr = 5177;

    // Receive halo from left rank
    MPI_Irecv(&fluid_particles[indexToReceiveLeft], num_from_left, Particletype, proc_to_left,tagl, MPI_COMM_WORLD, &edges->reqs[0]);
    // Receive halo from right rank
    MPI_Irecv(&fluid_particles[indexToReceiveRight], num_from_right, Particletype, proc_to_right,tagr, MPI_COMM_WORLD, &edges->reqs[1]);
    // Send halo to right rank
    MPI_Isend(sendr_buffer, num_moving_right, Particletype, proc_to_right, tagl, MPI_COMM_WORLD, &edges->reqs[2]);
    // Send halo to left rank
    MPI_Isend(sendl_buffer, num_moving_left, Particletype, proc_to_left, tagr, MPI_COMM_WORLD, &edges->reqs[3]);

    // Free allocated arays
    // Can't free send_buffer!!!!
    //!!!!
    // MEMORY LEAK HERE WITHOUT FREE, NEED TO FIX
    printf("MEMORY LEAK IN START HALO\n\n");
//    free(send_buffer);
}

void finishHaloExchange(fluid_particle *fluid_particles,  edge *edges, param *params)
{
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

    printf("rank %d, halo: recv %d from left, %d from right\n", params->rank,num_received_left,num_received_right);
}

// Transfer particles that are out of node bounds
void transferOOBParticles(fluid_particle *fluid_particles, oob *out_of_bounds, param *params)
{
    int i;
    fluid_particle *p;
    int rank = params->rank;
    int nprocs = params->nprocs;

    int i_left = 0;
    int i_right = 0;

    // Allocate send buffers
    printf("HACK MAX SEND NUMBER FOR NOW\n");
    int max_send = 1000;
    fluid_particle *send_buffer = malloc(max_send*sizeof(fluid_particle));
    fluid_particle *sendl_buffer = send_buffer;
    fluid_particle *sendr_buffer = send_buffer + max_send/2;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = &fluid_particles[i];

        // Set OOB particle indices and update number
        // remove particle from end of array to fill vacancy
        if (p->x < params->node_start_x && params->rank != 0) {
            sendl_buffer[i_left++] = *p;
            fluid_particles[i] = fluid_particles[params->number_fluid_particles_local-1];
            fluid_particles[i].id = i;
            params->number_fluid_particles_local--;
        }
        else if (p->x > params->node_end_x && params->rank != params->nprocs-1) {
            sendr_buffer[i_right++] = *p;
            fluid_particles[i] = fluid_particles[params->number_fluid_particles_local-1];
            fluid_particles[i].id = i;
            params->number_fluid_particles_local--;
        }
    }

    int num_moving_left = i_left;
    int num_moving_right = i_right;
    int total_send = num_moving_left + num_moving_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    // Get number of particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;
    int tag = 7006;
    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8278;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    // Allocate recieve buffers
    int total_recv = num_from_left + num_from_right;
    fluid_particle *recvl_buffer = &fluid_particles[params->number_fluid_particles_local];
    fluid_particle *recvr_buffer = &fluid_particles[params->number_fluid_particles_local + num_from_left];

    MPI_Status status;
    MPI_Request request;

    // Send oob particles to right processor receive oob particles from right processor
    int num_received_left = 0;
    int num_received_right = 0;

    // Sending to right, recv from left
    tag = 2522;
    MPI_Sendrecv(sendr_buffer, num_moving_right, Particletype,proc_to_right, tag,
                 recvl_buffer, num_from_left, Particletype, proc_to_left, tag,
                 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, Particletype, &num_received_left);

    // Sending to left, recv from right
    tag = 1165;
    MPI_Sendrecv(sendl_buffer, num_moving_left, Particletype, proc_to_left, tag,
                 recvr_buffer, num_from_right, Particletype, proc_to_right, tag,
                 MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, Particletype, &num_received_right);

    int total_received = num_received_right + num_received_left;

    // Update received particle ID
    for(i=params->number_fluid_particles_local; i<params->number_fluid_particles_local+total_received; i++)
    {
        fluid_particles[i].id = i;
    }

    // Update number of particles
    params->number_fluid_particles_local += total_received;

    printf("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n", rank, num_moving_left, num_moving_right, num_received_left, num_received_right);

    free(send_buffer);
}
