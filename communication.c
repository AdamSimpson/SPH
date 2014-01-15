#include <stdio.h>
#include "mpi.h"
#include "communication.h"
#include "fluid.h"
#include <stddef.h>

void createMpiTypes()
{
    //Create fluid particle type
    MPI_Datatype types[14];
    int i;
    for (i=0; i<12; i++) types[i] = MPI_DOUBLE;
    types[12] = MPI_INT;
    types[13] = MPI_UB;
    int blocklens[14];
    for (i=0; i<14; i++) blocklens[i] = 1;
    MPI_Aint disps[14];
    // Get displacement of each struct member
    disps[0] = offsetof( fluid_particle, x_prev);
    disps[1] = offsetof( fluid_particle, y_prev);
    disps[2] = offsetof( fluid_particle, x);
    disps[3] = offsetof( fluid_particle, y);
    disps[4] = offsetof( fluid_particle, v_x);
    disps[5] = offsetof( fluid_particle, v_y);
    disps[6] = offsetof( fluid_particle, a_x);
    disps[7] = offsetof( fluid_particle, a_y);
    disps[8] = offsetof( fluid_particle, density);
    disps[9] = offsetof( fluid_particle, density_near);
    disps[10] = offsetof( fluid_particle, pressure);
    disps[11] = offsetof( fluid_particle, pressure_near);
    disps[12] = offsetof( fluid_particle, id);
    disps[13] = sizeof(fluid_particle);
    // Commit type
    MPI_Type_create_struct( 14, blocklens, disps, types, &Particletype );
    MPI_Type_commit( &Particletype );
}

void freeMpiTypes()
{
    MPI_Type_free(&Particletype);
}

void startHaloExchange(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,  edge *edges, param *params)
{
    int i;
    int rank = params->rank;
    int nprocs = params->nprocs;

    fluid_particle *p;
    double h = params->smoothing_radius;

    // Set edge particle indicies and update number
    edges->number_edge_particles_left = 0;
    edges->number_edge_particles_right = 0;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = fluid_particle_pointers[i];
        if (p->x - params->node_start_x <= h)
            edges->edge_pointers_left[edges->number_edge_particles_left++] = p;
        else if (params->node_end_x - p->x <= h)
            edges->edge_pointers_right[edges->number_edge_particles_right++] = p;
    }

    int num_moving_left = edges->number_edge_particles_left;
    int num_moving_right = edges->number_edge_particles_right;

    // Setup nodes to left and right of self
    int proc_to_left =  (rank == 0 ? MPI_PROC_NULL : rank-1);
    int proc_to_right = (rank == nprocs-1 ? MPI_PROC_NULL : rank+1);

    printf("rank %d, halo: will send %d to left, %d to right\n", rank, num_moving_left, num_moving_right);

    // Get number of halo particles from right and left
    int num_from_left = 0;
    int num_from_right = 0;
    int tag = 3217;

    // Send number to right and receive from left
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8425;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("rank %d, halo: will recv %d from left, %d from right\n", rank, num_from_left, num_from_right);


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
   
    printf("rank %d, prams->max_fluid_particle_index: %d\n", rank,  params->max_fluid_particle_index);

    //Index to start receiving halo particles
    int indexToReceiveLeft = params->max_fluid_particle_index + 1;
    int indexToReceiveRight = indexToReceiveLeft + num_from_left;

    printf("rank %d, halo: send %d to left, %d to right, indexToRecieveLeft %d, indexToReceiveRight %d \n", rank, num_moving_left, num_moving_right, indexToReceiveLeft, indexToReceiveRight);

    int tagl = 4312;
    int tagr = 5177;
    // Receive halo from left rank
    MPI_Irecv(&fluid_particles[indexToReceiveLeft], num_from_left, Particletype, proc_to_left,tagl, MPI_COMM_WORLD, &edges->reqs[0]);
    // Receive halo from right rank
    MPI_Irecv(&fluid_particles[indexToReceiveRight], num_from_right, Particletype, proc_to_right,tagr, MPI_COMM_WORLD, &edges->reqs[1]);
    // Send halo to right rank
    MPI_Isend(fluid_particles,1,RightEdgetype,proc_to_right,tagl,MPI_COMM_WORLD, &edges->reqs[2]);
    MPI_Isend(fluid_particles,1,LeftEdgetype,proc_to_left,tagr,MPI_COMM_WORLD, &edges->reqs[3]);

    // Free allocated arays
    free(blocklens_left);
    free(blocklens_right);
    free(indicies_left);
    free(indicies_right);
}

void finishHaloExchange(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,  edge *edges, param *params)
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
void transferOOBParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, oob *out_of_bounds, param *params)
{
    int i;
    fluid_particle *p;
    int rank = params->rank;
    int nprocs = params->nprocs;

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
    MPI_Sendrecv(&num_moving_right, 1, MPI_INT, proc_to_right, tag, &num_from_left,1,MPI_INT,proc_to_left,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    // Send number to left and receive from right
    tag = 8278;
    MPI_Sendrecv(&num_moving_left, 1, MPI_INT, proc_to_left, tag, &num_from_right,1,MPI_INT,proc_to_right,tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

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
    MPI_Sendrecv(fluid_particles,1,RightMovetype,proc_to_right,tag,fluid_particles,1,LeftRecvtype,proc_to_left,tag,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status, Particletype, &num_received_left);
    // Sending to left, recv from right
    tag = 1165;
    MPI_Sendrecv(fluid_particles,1,LeftMovetype,proc_to_left,tag,fluid_particles,1,RightRecvtype,proc_to_right,tag,MPI_COMM_WORLD,&status);
    MPI_Get_count(&status, Particletype, &num_received_right);

    int total_sent = num_moving_left + num_moving_right;
    int total_received = num_received_right + num_received_left;

    printf("rank %d OOB: sent left %d, right: %d recv left:%d, right: %d\n", rank, num_moving_left, num_moving_right, num_received_left, num_received_right);
   
    // Update maximum particle index if neccessary
    int max_received_index = total_received?total_received-1:0;// If non received don't access indicies_recv[-1]...
    if (total_received && indicies_recv[max_received_index] > params->max_fluid_particle_index)
        params->max_fluid_particle_index = indicies_recv[max_received_index];
    
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

