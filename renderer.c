#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "circles_gl.h"
#include "mpi.h"
#include "communication.h"
#include "fluid.h"
#include "font_gl.h"

void start_renderer()
{
    // Setup initial OpenGL ES state
    // OpenGL state
    GL_STATE_T gl_state;
    memset(&gl_state, 0, sizeof(GL_STATE_T));

    // Start OpenGL
    init_ogl(&gl_state);

    // Initialize circle OpenGL state
    CIRCLE_T circle_state;
    init_circles(&circle_state);

    // Initialize font atlas
    FONT_T font_state;
    init_font(&font_state, gl_state.screen_width, gl_state.screen_height);

    // enum to handle currently selected parameter
    parameters selected_param;
    // Set to first parameter
    selected_param = 0;

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

    printf("nprocs %d\n", params[0].nprocs);

    // Allocate particle receive array
    int max_particles = params[0].number_fluid_particles_global;
    int num_coords = 2;
    float *particle_coords = (float*)malloc(num_coords * max_particles*sizeof(float));

    // Allocate points array(position + color)
    int point_size = 5 * sizeof(float);
    float *points = (float*)malloc(point_size*max_particles);

    // Allocate mover point array(position + color)
    float mover_point[5];

    // Particle coordinate gather values
    int *particle_counts = malloc(num_procs * sizeof(int));
    int *particle_displs = malloc(num_procs * sizeof(int));

    // Set background color
    glClearColor(0.121f, 0.243f, 0.423f, 1.0f);

    // Create color index - hard coded for now to experiment
    float colors_by_rank[9] = {0.69,0.07,0.07,
        1.0,1.0,0.1,
        0.08,0.52,0.8};

    // Perhaps the RECV loop will help pipeline particle send and draw more than a gather
    int num_coords_rank;
    int total_coords, current_rank;
    double mouse_x, mouse_y, mouse_x_scaled, mouse_y_scaled;
    double mover_radius, mover_radius_scaled;

    int frames_per_fps = 30;
    int num_steps = 0;
    double current_time;
    double wall_time = MPI_Wtime();
    double fps;

    while(1){
	// Every frames_per_fps steps calculate FPS
	if(num_steps%frames_per_fps == 0) {
	    current_time =  MPI_Wtime();
	    wall_time = current_time - wall_time;
	    fps = frames_per_fps/wall_time;
	    num_steps = 0;
	    wall_time = current_time;
	}	    

	// Get keyboard key press
	// process appropriately
	check_key_press(&gl_state);	

        // Recieve paramaters struct from all nodes
        MPI_Gatherv(MPI_IN_PLACE, 0, Paramtype, params, param_counts, param_displs, Paramtype, 0, MPI_COMM_WORLD);

        // Update paramaters as needed
        get_mouse(&mouse_x, &mouse_y, &gl_state);
        mouse_x_scaled = mouse_x*10.0 + 10.0;
        mouse_y_scaled = mouse_y*10.0 + 10.0;
        mover_radius = 1.0;
        for(i=0; i< num_compute_procs; i++) {
            params[i].mover_center_x = mouse_x_scaled;
            params[i].mover_center_y = mouse_y_scaled;
            params[i].mover_radius = mover_radius;
        }

        // Send updated paramaters to compute nodes
        MPI_Scatterv(params, param_counts, param_displs, Paramtype, MPI_IN_PLACE, 0, Paramtype, 0, MPI_COMM_WORLD);

        // Retrieve all particle coordinates (x,y)
        num_coords_rank = 0;
        total_coords = 0;
        for(i=0; i<num_procs; i++) {
            num_coords_rank = i?2*params[i-1].number_fluid_particles_local:0;
            particle_counts[i] = i?num_coords_rank:0;
            particle_displs[i] = i?total_coords:0;
            total_coords += num_coords_rank;
        }
        MPI_Gatherv(MPI_IN_PLACE, 0, MPI_FLOAT, particle_coords, particle_counts, particle_displs, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Create points array (x,y,r,g,b)
        current_rank = 0;
        int particle_count = 1;
        for(j=0; j<total_coords/2; j++, particle_count+=2) {
            if ( particle_count > particle_counts[1+current_rank]){
                current_rank++;
                particle_count = 1;
            }

            points[j*5]   = particle_coords[j*2]/10.0 - 1.0; 
            points[j*5+1] = particle_coords[j*2+1]/10.0 - 1.0;
            points[j*5+2] = colors_by_rank[3*current_rank];
            points[j*5+3] = colors_by_rank[3*current_rank+1];
            points[j*5+4] = colors_by_rank[3*current_rank+2];
        }

        // Clear background
        glClear(GL_COLOR_BUFFER_BIT);

        // Render particles
        update_points(points, total_coords/2, &circle_state);

        // Render mover
        mover_point[0] = mouse_x_scaled/10.0 - 1.0;
        mover_point[1] = mouse_y_scaled/10.0 - 1.0;
        mover_point[2] = 0.0;
        mover_point[3] = 0.0;
        mover_point[4] = 0.0;
        mover_radius_scaled = mover_radius*70.0;
        update_mover_point(mover_point, mover_radius_scaled, &circle_state);

	    // Draw font parameters
//        render_parameters(&font_state, selected_param, params[0].g, 1.0, 1.0, 1.0, 1.0);

        // Draw FPS
  	    render_fps(&font_state, fps);

        // Swap front/back buffers
        swap_ogl(&gl_state);

        num_steps++;
    }

//    exit_ogl(&state.gl_state);

}
