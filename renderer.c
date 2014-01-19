#include "circles_gl.h"
#include "stdio.h"
#include "mpi.h"
#include <string.h> 
#include <stdlib.h>

void start_renderer()
{
    // Setup initial OpenGL ES state
    STATE_T state;
    memset(&state, 0, sizeof(STATE_T));

    // Start OGLES
    init_ogl(&state.gl_state);

    // Create and set shaders
    create_shaders(&state);

    // Allocate particle receive array
    int max_particles = 5000;
    int num_coords = 2;
    float *positions = (float*)malloc(num_coords * max_particles*sizeof(float));

    // Allocate points array(position + color);
    int point_size = 5 * sizeof(float);
    float *points = (float*)malloc(point_size*max_particles);

    int i,j, coords_recvd, disp;
    int num_compute = 3;
    MPI_Status status;

    // Perhaps the RECV loop will help pipeline particle send and draw more than a gather
    while(1){

        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

	for(i=0; i<num_compute; i++) {
	    // receive particles
	    MPI_Recv(positions, max_particles, MPI_FLOAT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_FLOAT, &coords_recvd);
	    coords_recvd/=2;

            // Create proper points array
            for(j=0; j<coords_recvd; j++) {
	        points[j*5]   = positions[j*2]/10.0 - 1.0; 
	        points[j*5+1] = positions[j*2+1]/10.0 - 0.8;
                points[j*5+2] = 0.0;
	        points[j*5+3] = 0.0;
		points[j*5+4] = 1.0;
	
            }

  	    // Render particles
            update_points(points, coords_recvd, &state);
        }

	// Swap front/back buffers
        swap_ogl(&state.gl_state);
    };

    exit_ogl(&state.gl_state);

}
