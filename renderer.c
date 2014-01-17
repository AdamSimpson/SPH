#include "egl_utils.h"
#include "circles_gl.h"
#include "stdio.h"
#include "linux/input.h"
#include "mpi.h"

void start_renderer()
{
    // Setup initial OpenGL ES state
    STATE_T state;
    memset(&state, 0, sizeof(STATE_T));

    bcm_host_init();

    // Start OGLES
    init_ogl(&state.egl_state);

    // Create buffer
    glGenBuffers(1, &state.vbo);

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
    int num_compute = 2;
    MPI_Status status;

    // Perhaps the RECV loop will help pipeline particle send and draw more than a gather
    while(1){
	for(i=0; i<num_compute; i++) {
	    // receive particles
	    MPI_Recv(positions, max_particles, MPI_FLOAT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_FLOAT, &coords_recvd);
	    coords_recvd/=2;

            // Create proper points array
            for(j=0; j<coords_recvd; j++) {
	        points[j*5]   = positions[j*2]/20.0 - 0.5; 
	        points[j*5+1] = positions[j*2+1]/10.0 - 0.5;
                points[j*5+2] = 0.0;
	        points[j*5+3] = 0.0;
		points[j*5+4] = 1.0;
            }

  	    // Render particles
            update_points(points, coords_recvd, &state);
        }

	// Swap buffers
        egl_swap(&state.egl_state);
    };

    exit_ogl(&state.egl_state);

}
