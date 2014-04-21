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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "particles_gl.h"
#include "mover_gl.h"
#include "background_gl.h"
#include "mpi.h"
#include "geometry.h"
#include "communication.h"
#include "fluid.h"
#include "font_gl.h"
#include "dividers_gl.h"
#include "renderer.h"

#ifdef LIGHT
    #include "rgb_light.h"
#endif

void start_renderer()
{
    // Setup initial OpenGL state
    gl_t gl_state;
    memset(&gl_state, 0, sizeof(gl_t));

    // Start OpenGL
    render_t render_state;
    init_ogl(&gl_state, &render_state);
    render_state.show_dividers = false;
    render_state.pause = false;
    render_state.screen_width = gl_state.screen_width;
    render_state.screen_height = gl_state.screen_height;

    // Initialize particles OpenGL state
    particles_t particle_GLstate;
    init_particles(&particle_GLstate, gl_state.screen_width, gl_state.screen_height);

    // Initialize mover OpenGL state
    mover_t mover_GLstate;
    init_mover(&mover_GLstate);

    // Initialize font OpenGL state
    font_t font_state;
    init_font(&font_state, gl_state.screen_width, gl_state.screen_height);

    // Initialize background OpenGL state
    background_t background_state;
    init_background(&background_state, gl_state.screen_width, gl_state.screen_height);

    // Initialize node divider OpenGL state
    dividers_t dividers_state;
    init_dividers(&dividers_state, gl_state.screen_width, gl_state.screen_height);

    // Initialize RGB Light if present
    #ifdef LIGHT
    rgb_light_t light_state;
    init_rgb_light(&light_state, 255, 0, 0);
    #endif

    // Number of processes
    int num_procs, num_compute_procs, num_compute_procs_active;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    num_compute_procs = num_procs - 1;
    num_compute_procs_active = num_compute_procs;

    // Allocate array of paramaters
    // So we can use MPI_Gather instead of MPI_Gatherv
    tunable_parameters *node_params = malloc(num_compute_procs*sizeof(tunable_parameters));

    // The render node must keep it's own set of master parameters
    // This is due to the GLFW key callback method
    tunable_parameters *master_params = malloc(num_compute_procs*sizeof(tunable_parameters));

    // Setup render state
    render_state.node_params = node_params;
    render_state.master_params = master_params;
    render_state.num_compute_procs = num_compute_procs;
    render_state.num_compute_procs_active = num_compute_procs;
    render_state.selected_parameter = 0;

    int i,j;

    // Broadcast pixels ratio
    short pixel_dims[2];
    pixel_dims[0] = (short)gl_state.screen_width;
    pixel_dims[1] = (short)gl_state.screen_height;
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
 
    // Recv simulation world dimensions from global rank 1
    float sim_dims[2];
    MPI_Recv(sim_dims, 2, MPI_FLOAT, 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    render_state.sim_width = sim_dims[0];
    render_state.sim_height = sim_dims[1];
    // Receive number of global particles
    int max_particles;
    MPI_Recv(&max_particles, 1, MPI_INT, 1, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Calculate world unit to pixel
    float world_to_pix_scale = gl_state.screen_width/render_state.sim_width;

    // Gatherv initial tunable parameters values
    int *param_counts = malloc(num_procs * sizeof(int));
    int *param_displs = malloc(num_procs * sizeof(int));
    for(i=0; i<num_procs; i++) {
        param_counts[i] = i?1:0; // will not receive from rank 0
        param_displs[i] = i?i-1:0; // rank i will reside in params[i-1]
    }
    // Initial gather
    MPI_Gatherv(MPI_IN_PLACE, 0, TunableParamtype, node_params, param_counts, param_displs, TunableParamtype, 0, MPI_COMM_WORLD);

    // Fill in master parameters
    for(i=0; i<render_state.num_compute_procs; i++)
        render_state.master_params[i] = node_params[i];

    // Set mover state
    mover_GLstate.mover_type = render_state.master_params[0].mover_type;

    // Allocate particle receive array
    int num_coords = 2;
    short *particle_coords = malloc(num_coords * max_particles*sizeof(short));

    // Allocate points array(position + color)
    int point_size = 5 * sizeof(float);
    float *points = malloc(point_size*max_particles);

    // Allocate mover point array(position + color)
    float mover_center[2];
    float mover_color[3];

    // Allocate space for node edges
    float *node_edges = (float*)malloc(2*render_state.num_compute_procs*sizeof(float));

    // Number of coordinates received from each proc
    int *particle_coordinate_counts = malloc(num_compute_procs * sizeof(int));
    // Keep track of order in which particles received
    int *particle_coordinate_ranks = malloc(num_compute_procs * sizeof(int));

    // Set background color
    glClearColor(0.15, 0.15, 0.15, 1.0);

    // Create color index, equally spaced around HSV
    float *colors_by_rank = malloc(3*render_state.num_compute_procs*sizeof(float));
    float angle_space = 0.5f/(float)render_state.num_compute_procs;
    float HSV[3];
    for(i=0; i<render_state.num_compute_procs; i++)
    {
        if(i%2)
            HSV[0] = angle_space*i;
        else
            HSV[0] = angle_space*i + 0.5f;
        HSV[1] = 1.0f;
        HSV[2] = 0.8f;
        hsv_to_rgb(HSV, colors_by_rank+3*i);
    }
 
    #ifdef LIGHT
    MPI_Bcast(colors_by_rank, 3*render_state.num_compute_procs, MPI_FLOAT, 0, MPI_COMM_WORLD);
    #endif

    int num_coords_rank;
    int current_rank, num_parts;
    float mover_gl_dims[2];

    int frames_per_fps = 30;
    #ifdef RASPI
    int frames_per_check = 1;
    #else
    int frames_per_check = 1;
    #endif
    int num_steps = 0;
    double current_time;
    double wall_time = MPI_Wtime();
    float fps=0.0f;

    // Setup MPI requests used to gather particle coordinates
    MPI_Request coord_reqs[num_compute_procs];
    int src, coords_recvd;
    float gl_x, gl_y;
    // Particle radius in pixels
    float particle_diameter_pixels = gl_state.screen_width * 0.0125;

    MPI_Status status;

    // Remove all partitions but one initially
    for(i=0; i<render_state.num_compute_procs-1; i++)
        remove_partition(&render_state);

    while(1){
        // Every frames_per_fps steps calculate FPS
        if(num_steps%frames_per_fps == 0) {
            current_time =  MPI_Wtime();
            wall_time = current_time - wall_time;
            fps = frames_per_fps/wall_time;
            num_steps = 0;
            wall_time = current_time;
        }

        // Check to see if simulation should close
        if(window_should_close(&gl_state)) {
            for(i=0; i<render_state.num_compute_procs; i++)
                render_state.node_params[i].kill_sim = true;
            // Send kill paramaters to compute nodes
            MPI_Scatterv(node_params, param_counts, param_displs, TunableParamtype, MPI_IN_PLACE, 0, TunableParamtype, 0, MPI_COMM_WORLD);
            break;
        }    

        // Check for user keyboard/mouse input
        if(render_state.pause) {
            while(render_state.pause)
                check_user_input(&gl_state);
        }
        else
            check_user_input(&gl_state);

        // Update node params with master param values
        update_node_params(&render_state);

        // Send updated paramaters to compute nodes
        MPI_Scatterv(node_params, param_counts, param_displs, TunableParamtype, MPI_IN_PLACE, 0, TunableParamtype, 0, MPI_COMM_WORLD);

        // Retrieve all particle coordinates (x,y)
  	    // Potentially probe is expensive? Could just allocated num_compute_procs*num_particles_global and async recv
	    // OR do synchronous recv...very likely that synchronous receive is as fast as anything else
	    coords_recvd = 0;
	    for(i=0; i<render_state.num_compute_procs; i++) {
	        // Wait until message is ready from any proc
                MPI_Probe(MPI_ANY_SOURCE, 17, MPI_COMM_WORLD, &status);
	        // Retrieve probed values
                src = status.MPI_SOURCE;
                particle_coordinate_ranks[i] = src-1;
	        MPI_Get_count(&status, MPI_SHORT, &particle_coordinate_counts[src-1]); // src-1 to account for render node
	        // Start async recv using probed values
	        MPI_Irecv(particle_coords + coords_recvd, particle_coordinate_counts[src-1], MPI_SHORT, src, 17, MPI_COMM_WORLD, &coord_reqs[src-1]);
            // Update total number of floats recvd
            coords_recvd += particle_coordinate_counts[src-1];
	    }

        // Ensure a balanced partition
        // We pass in number of coordinates instead of particle counts    
        if(num_steps%frames_per_check == 0)
            check_partition_left(&render_state, particle_coordinate_counts, coords_recvd);

        // Clear background
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw background image
        draw_background(&background_state);

        // update mover
        sim_to_opengl(&render_state, render_state.master_params[0].mover_center_x, render_state.master_params[0].mover_center_y, &gl_x, &gl_y);
        mover_center[0] = gl_x;
        mover_center[1] = gl_y;
        mover_color[0] = 1.0f;
        mover_color[1] = 0.0f;
        mover_color[2] = 0.0f;
        // Mover bounding rectangle half width/height lengths in ogl system
        // Subtract off particle diamter so no particle/mover penetration
        mover_gl_dims[0] = render_state.master_params[0].mover_width/(render_state.sim_width*0.5f) - particle_diameter_pixels/(gl_state.screen_width*0.5f) ;
        mover_gl_dims[1] = render_state.master_params[0].mover_height/(render_state.sim_height*0.5f) - particle_diameter_pixels/(gl_state.screen_height*0.5f);

        render_all_text(&font_state, &render_state, fps);

        if(render_state.show_dividers)
        {
            // Convert node start/end parameters into GL coordinates and render
            for(i=0; i<render_state.num_compute_procs_active; i++)
            {
                float start_gl_x, end_gl_x;
                float null_y;
                sim_to_opengl(&render_state, node_params[i].node_start_x, 0.0, &start_gl_x, &null_y);
                sim_to_opengl(&render_state, node_params[i].node_end_x, 0.0, &end_gl_x, &null_y);
                node_edges[2*i] = start_gl_x;
                node_edges[2*i+1] = end_gl_x;
            }
            render_dividers(&dividers_state, node_edges, colors_by_rank, render_state.num_compute_procs_active);
        }

        // Wait for all coordinates to be received
        MPI_Waitall(num_compute_procs, coord_reqs, MPI_STATUSES_IGNORE);

        // Create points array (x,y,r,g,b)
	    i = 0;
        current_rank = particle_coordinate_ranks[i];
        // j == coordinate pair
        for(j=0, num_parts=1; j<coords_recvd/2; j++, num_parts++) {
	    // Check if we are processing a new rank's particles
            if ( num_parts > particle_coordinate_counts[current_rank]/2){
                current_rank =  particle_coordinate_ranks[++i];
                num_parts = 1;
		// Find next rank with particles if current_rank has 0 particles
		while(!particle_coordinate_counts[current_rank])
                    current_rank = particle_coordinate_ranks[++i];
            }

            points[j*5]   = particle_coords[j*2]/(float)SHRT_MAX; 
            points[j*5+1] = particle_coords[j*2+1]/(float)SHRT_MAX;
            points[j*5+2] = colors_by_rank[3*current_rank];
            points[j*5+3] = colors_by_rank[3*current_rank+1];
            points[j*5+4] = colors_by_rank[3*current_rank+2];

        }

	    // Draw particles
        render_particles(points, particle_diameter_pixels, coords_recvd/2, &particle_GLstate);

        // Render over particles to hide penetration
        render_mover(mover_center, mover_gl_dims, mover_color, &mover_GLstate);

        // Swap front/back buffers
        swap_ogl(&gl_state);

        num_steps++;
    }

    // Clean up memory
    exit_ogl(&gl_state);
    free(node_params);
    free(master_params);
    free(param_counts);
    free(param_displs);
    free(particle_coords);
    free(points);
    free(particle_coordinate_counts);
    free(particle_coordinate_ranks);
    free(colors_by_rank);
}

// Translate between OpenGL coordinates with origin at screen center
// to simulation coordinates
void opengl_to_sim(render_t *render_state, float x, float y, float *sim_x, float *sim_y)
{
    float half_width = render_state->sim_width*0.5f;
    float half_height = render_state->sim_height*0.5f;

    *sim_x = x*half_width + half_width;
    *sim_y = y*half_height + half_height;
}

// Translate between simulation coordinates, origin bottom left, and opengl -1,1 center of screen coordinates
void sim_to_opengl(render_t *render_state, float x, float y, float *gl_x, float *gl_y)
{
    float half_width = render_state->sim_width*0.5f;
    float half_height = render_state->sim_height*0.5f;

    *gl_x = x/half_width - 1.0f;
    *gl_y = y/half_height - 1.0f;
}

void update_node_params(render_t *render_state)
{
    int i;
	// Update all node parameters with master paramter values
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->node_params[i] = render_state->master_params[i]; 
}

// Checks for a balanced number of particles on each compute node
// If unbalanced the partition between nodes will change 
// Check from right to left
void check_partition_left(render_t *render_state, int *particle_counts, int total_particles)
{
    int rank, diff;
    float h, dx, length, length_left, length_right;   

    // Particles per proc if evenly divided
    int even_particles = total_particles/render_state->num_compute_procs_active;
    int max_diff = even_particles/15.0f;

    // Fixed distance to move partition is 0.125*smoothing radius
    h = render_state->master_params[0].smoothing_radius;
    dx = h*0.125;

    tunable_parameters *master_params = render_state->master_params;

    for(rank=render_state->num_compute_procs_active; rank-- > 1; )
    {
        length =  master_params[rank].node_end_x - master_params[rank].node_start_x;
        length_left =  master_params[rank-1].node_end_x - master_params[rank-1].node_start_x;
        diff = particle_counts[rank] - even_particles;

        // current rank has too many particles
        if( diff > max_diff && length > 2*h) {
            master_params[rank].node_start_x += dx;
            master_params[rank-1].node_end_x = master_params[rank].node_start_x;
        }
        // current rank has too few particles
        else if (diff < -max_diff && length_left > 2*h) {
            master_params[rank].node_start_x -= dx;
            master_params[rank-1].node_end_x = master_params[rank].node_start_x;
        }
    }

    // Left most partition has a tendency to not balance correctly so we test it explicitly
    length =  master_params[0].node_end_x - master_params[0].node_start_x;
    length_right =  master_params[1].node_end_x - master_params[1].node_start_x;
    diff = particle_counts[0] - even_particles;

    // current rank has too many particles
    if( diff > max_diff && length > 2*h) {
        master_params[0].node_end_x -= dx;
        master_params[1].node_start_x = master_params[0].node_end_x;
    }
    else if (diff < -max_diff && length_right > 2*h) {
        master_params[0].node_end_x += dx;
        master_params[1].node_start_x = master_params[0].node_end_x;
    }
}

// Convert hsv to rgb
// input hsv [0:1]
// output rgb [0,1]
void hsv_to_rgb(float* HSV, float *RGB)
{
    float hue, saturation, value, hh, p, q, t, ff, r, g, b;
    long i;

    hue = HSV[0];
    saturation = HSV[1];
    value = HSV[2];

    hh = hue*360.0f;
    if(hh >= 360.0f)
	hh = 0.0f;
    hh /= 60.0f;
    i = (long)hh;
    ff = hh - i;
    p = value * (1.0f - saturation);
    q = value * (1.0f - (saturation * ff));
    t = value * (1.0f - (saturation * (1.0f - ff)));

    switch(i) {
        case 0:
	    r = value;
	    g = t;
	    b = p;
	    break;
	case 1:
	    r = q;
	    g = value;
	    b = p;
	    break;
	case 2:
	    r = p;
	    g = value;
	    b = t;
	    break;
	case 3:
	    r = p;
	    g = q;
	    b = value;
	    break;
        case 4:
	    r = t;
	    g = p;
	    b = value;
	    break;
	case 5:
	    r = value;
	    g = p;
	    b = q;
	    break;
	default:
	    r = value;
	    g = p;
	    b = q;
	    break;
    }

    RGB[0] = r;
    RGB[1] = g;
    RGB[2] = b;
}
