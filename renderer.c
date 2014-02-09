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
    RENDER_T render_state;
    init_ogl(&gl_state, &render_state);

    // Initialize circle OpenGL state
    CIRCLE_T circle_state;
    init_circles(&circle_state);

    // Initialize font atlas
    FONT_T font_state;
    init_font(&font_state, gl_state.screen_width, gl_state.screen_height);

    // Number of processes
    int num_procs, num_compute_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    num_compute_procs = num_procs - 1;

    // Allocate array of paramaters
    // So we can use MPI_Gather instead of MPI_Gatherv
    param *params = malloc(num_compute_procs*sizeof(param));

    // Setup render state
    render_state.params = params;
    render_state.num_params = num_compute_procs;
    render_state.selected_parameter = 0;

    int i,j;

    // Broadcast aspect ratio
    float aspect_ratio = (float)gl_state.screen_width/(float)gl_state.screen_height;
    MPI_Bcast(&aspect_ratio, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
 
    // Recv world dimensions from global rank 1
    float world_dims[2];
    MPI_Recv(world_dims, 2, MPI_FLOAT, 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Calculate world unit to pixel
    float world_to_pix_scale = gl_state.screen_width/world_dims[0];

    // Gatherv values
    int *param_counts = malloc(num_procs * sizeof(int));
    int *param_displs = malloc(num_procs * sizeof(int));
    for(i=0; i<num_procs; i++) {
        param_counts[i] = i?1:0; // will not receive from rank 0
        param_displs[i] = i?i-1:0; // rank i will reside in params[i-1]
    }
    // Initial gather
    MPI_Gatherv(MPI_IN_PLACE, 0, Paramtype, params, param_counts, param_displs, Paramtype, 0, MPI_COMM_WORLD);

    printf("global particles: %d\n", params[0].number_fluid_particles_global);

    // Allocate particle receive array
    int max_particles = params[0].number_fluid_particles_global;
    int num_coords = 2;
    float *particle_coords = (float*)malloc(num_coords * max_particles*sizeof(float));

    // Allocate points array(position + color)
    int point_size = 5 * sizeof(float);
    float *points = (float*)malloc(point_size*max_particles);

    // Allocate mover point array(position + color)
    float mover_point[5];

    // Number of coordinates received from each proc
    int *particle_coordinate_counts = malloc(num_compute_procs * sizeof(int));

    // Set background color
    glClearColor(0, 0, 0, 1);

    // Create color index - hard coded for now to experiment
    float colors_by_rank[9] = {0.69,0.07,0.07,
        1.0,1.0,0.1,
        0.08,0.52,0.8};

    int num_coords_rank;
    int current_rank;
    float mouse_x, mouse_y, mouse_x_scaled, mouse_y_scaled;
    float mover_radius, mover_radius_scaled;

    int frames_per_fps = 30;
    int num_steps = 0;
    double current_time;
    double wall_time = MPI_Wtime();
    float fps=0.0f;

    // Setup MPI requests used to gather particle coordinates
    MPI_Request coord_reqs[num_compute_procs];
    MPI_Status status;
    int src, coords_recvd;
    float gl_x, gl_y;
    float particle_radius = 1.0f;

    while(1){
	// Every frames_per_fps steps calculate FPS
	if(num_steps%frames_per_fps == 0) {
	    current_time =  MPI_Wtime();
	    wall_time = current_time - wall_time;
	    fps = frames_per_fps/wall_time;
	    num_steps = 0;
	    wall_time = current_time;
	}	    

        // Send updated paramaters to compute nodes
//        MPI_Scatterv(params, param_counts, param_displs, Paramtype, MPI_IN_PLACE, 0, Paramtype, 0, MPI_COMM_WORLD);

        // Get keyboard key press
        // process appropriately
        check_key_press(&gl_state);

        // Update mover position
        get_mouse(&mouse_x, &mouse_y, &gl_state);
        pixel_to_sim(world_dims, mouse_x, mouse_y, &mouse_x_scaled, &mouse_y_scaled);
        mover_radius = 1.0f;
        render_state.mover_center_x = mouse_x_scaled;
        render_state.mover_center_y = mouse_y_scaled;
        render_state.mover_radius = mover_radius;

        // Retrieve all particle coordinates (x,y)
	// Potentially probe is expensive? Could just allocated num_compute_procs*num_particles_global and async recv
	// OR do synchronous recv
	coords_recvd = 0;
	for(i=0; i<num_compute_procs; i++) {
	    // Wait until message is ready from any proc
            MPI_Probe(MPI_ANY_SOURCE, 17, MPI_COMM_WORLD, &status);
	    // Retrieve probed values
    	    src = status.MPI_SOURCE;
	    MPI_Get_count(&status, MPI_FLOAT, &particle_coordinate_counts[src-1]); // src-1 to account for render node
	    // Start async recv using probed values
	    MPI_Irecv(particle_coords + coords_recvd, particle_coordinate_counts[src-1], MPI_FLOAT, src, 17, MPI_COMM_WORLD, &coord_reqs[src-1]);
            // Update total number of floats recvd
            coords_recvd += particle_coordinate_counts[src-1];
	}

        // Clear background
        glClear(GL_COLOR_BUFFER_BIT);

        // Render mover
        sim_to_opengl(world_dims, mouse_x_scaled, mouse_y_scaled, &gl_x, &gl_y);
        mover_point[0] = gl_x;
        mover_point[1] = gl_y;
        mover_point[2] = 1.0f;
        mover_point[3] = 1.0f;
        mover_point[4] = 1.0f;
        mover_radius_scaled = mover_radius*world_to_pix_scale - particle_radius;
        update_mover_point(mover_point, mover_radius_scaled, &circle_state);

        // Draw FPS
        render_fps(&font_state, fps);

        // Draw font parameters
        // SHOULD JUST PASS IN render_state
        // RENDER STATE SHOULD INCLUDE PARAMETER VALUES TO DISPLAY
        render_parameters(&font_state, render_state.selected_parameter, render_state.g, 1.0f, 1.0f, 1.0f, 1.0f);

	// Wait for all coordinates to be received
	MPI_Waitall(num_compute_procs, coord_reqs, MPI_STATUSES_IGNORE);

/*
        // Create points array (x,y,r,g,b)
        current_rank = 0;
        int particle_count = 1;
        for(j=0; j<coords_recvd/2; j++, particle_count+=2) {
            if ( particle_count > particle_coordinate_counts[1+current_rank]){
                current_rank++;
                particle_count = 1;
            }
            sim_to_opengl(world_dims, particle_coords[j*2], particle_coords[j*2+1], &gl_x, &gl_y);
            points[j*5]   = gl_x; 
            points[j*5+1] = gl_y;
            points[j*5+2] = colors_by_rank[3*current_rank];
            points[j*5+3] = colors_by_rank[3*current_rank+1];
            points[j*5+4] = colors_by_rank[3*current_rank+2];
        }

	// Draw particles
        update_points(points, particle_radius, coords_recvd/2, &circle_state);
*/
        // Swap front/back buffers
        swap_ogl(&gl_state);

        num_steps++;

        // TODO: this function!
        // Calculate problem partitioning
//        balance_partitions(render_state);
    }

//    exit_ogl(&state.gl_state);

}

// Move selected parameter up
void move_parameter_up(RENDER_T *render_state)
{
    if(render_state->selected_parameter == MIN)
        render_state->selected_parameter = MAX;
    else
	render_state->selected_parameter--;
}

// Move selected parameter down
void move_parameter_down(RENDER_T *render_state) 
{
    if(render_state->selected_parameter == MAX)
        render_state->selected_parameter = MIN;
    else
        render_state->selected_parameter++;
}

void increase_parameter(RENDER_T *render_state)
{
    switch(render_state->selected_parameter) {
        case GRAVITY:
	   increase_gravity(render_state);
	    break;
    }
}

void decrease_parameter(RENDER_T *render_state)
{
    switch(render_state->selected_parameter) {
        case GRAVITY:
            decrease_gravity(render_state);
            break;
    }

}

// Increase gravity parameter
void increase_gravity(RENDER_T *render_state)
{
    static const float max_grav = -9.0;
    if(render_state->params[0].g < max_grav)
        return;

    int i;
    for(i=0; i<render_state->num_params; i++) {
        render_state->params[i].g -= 1.0;
    }
}

// Decreate gravity parameter
void decrease_gravity(RENDER_T *render_state)
{
    static const float min_grav = 9.0;
    if(render_state->params[0].g > min_grav)
        return;

    int i;
    for(i=0; i<render_state->num_params; i++)
        render_state->params[i].g += 1.0;
}

// Translate between pixel coordinates with origin at screen center
// to simulation coordinates
void pixel_to_sim(float *world_dims, float x, float y, float *sim_x, float *sim_y)
{
    float half_width = world_dims[0]*0.5;
    float half_height = world_dims[1]*0.5;

    *sim_x = x*half_width + half_width;
    *sim_y = y*half_height + half_height;
}

void sim_to_opengl(float *world_dims, float x, float y, float *gl_x, float *gl_y)
{
    float half_width = world_dims[0]*0.5;
    float half_height = world_dims[1]*0.5;

    *gl_x = x/half_width - 1.0;
    *gl_y = y/half_height - 1.0;
}
