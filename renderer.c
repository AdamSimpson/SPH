#include "circles_gl.h"
#include "stdio.h"
#include "mpi.h"
#include "communication.h"
#include "fluid.h"
#include <string.h> 
#include <stdlib.h>

void start_renderer()
{
    // Setup initial OpenGL ES state
    STATE_T state;
    memset(&state, 0, sizeof(STATE_T));

    // Start OpenGL
    init_ogl(&state.gl_state);

    // Create OpenGL buffers
    create_buffers(&state);

    // Create and set shaders
    create_shaders(&state);

    // Number of processes
    int num_procs, num_compute_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    num_compute_procs = num_procs - 1;

    // Allocate array of paramaters
    // So we can use MPI_Gather instead of MPI_Gatherv
    param *params = malloc(num_compute_procs*sizeof(param));

    int i,j;

    // Gatherv values
    int *param_counts = malloc(num_procs * sizeof(int));
    int *param_displs = malloc(num_procs * sizeof(int));
    for(i=0; i<num_procs; i++) {
        param_counts[i] = i?1:0; // will not receive from rank 0
        param_displs[i] = i?i-1:0; // rank i will reside in params[i-1]
    }

    // Initial gather
    MPI_Gatherv(MPI_IN_PLACE, 0, Paramtype, params, param_counts, param_displs, Paramtype, 0, MPI_COMM_WORLD);

    printf("global on render: %d\n", params[0].number_fluid_particles_global);

    // Allocate particle receive array
    int max_particles = params[0].number_fluid_particles_global;
    int num_coords = 2;
    float *particle_coords = (float*)malloc(num_coords * max_particles*sizeof(float));

    // Allocate points array(position + color);
    int point_size = 5 * sizeof(float);
    float *points = (float*)malloc(point_size*max_particles);

    int coords_recvd, disp;
    MPI_Status status;

    // Particle coordinate gather values
    int *particle_counts = malloc(num_procs * sizeof(int));
    int *particle_displs = malloc(num_procs * sizeof(int));

    // Perhaps the RECV loop will help pipeline particle send and draw more than a gather
    int num_coords_rank;
    int total_coords;
    while(1){

        // Recieve paramaters struct from all nodes
        MPI_Gatherv(MPI_IN_PLACE, 0, Paramtype, params, param_counts, param_displs, Paramtype, 0, MPI_COMM_WORLD);

        // Update paramaters as needed

        // Send updated paramaters to compute nodes
        MPI_Scatterv(params, param_counts, param_displs, Paramtype, MPI_IN_PLACE, 0, Paramtype, 0, MPI_COMM_WORLD);

        // Retrieve all particle coordinates (x,y)
        num_coords_rank = 0;
        total_coords = 0;
	for(i=0; i<num_procs; i++) {
	    num_coords_rank = i?2*params[i-1]->number_fluid_particles_local:0;
	    total_coords += i?num_coords_rank:0;
            particle_counts[i] = i?num_coords_rank:0;
            particle_displs[i] = i?total_coords:0;
        }
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_FLOAT, particle_coords, particle_counts, particle displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Create proper points array
        for(j=0; j<coords_recvd; j++) {
	    points[j*5]   = positions[j*2]/10.0 - 1.0; 
            points[j*5+1] = positions[j*2+1]/10.0 - 0.8;
            points[j*5+2] = 0.0;
	    points[j*5+3] = 0.0;
            points[j*5+4] = 1.0;
        }

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Render particles
        update_points(points, coords_recvd, &state);

	// Swap front/back buffers
        swap_ogl(&state.gl_state);
    }

    exit_ogl(&state.gl_state);

}
