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
#include "mpi.h"
#include "setup.h"
#include "mover_gl.h"
#include "communication.h"
#include "fluid.h"
#include "font_gl.h"
#include "container_gl.h"
#include "camera.hpp"
#include "renderer.hpp"

#ifdef BLINK1
    #include "blink1_light.h"
#endif

void Renderer::start_rendering()
{
    int i,j;

    this->set_activity_time();

    // Initialize particles OpenGL state
    particles_t particle_GLstate;
    init_particles(&particle_GLstate, this->screen_width(), this->screen_height());

    // Initialize mover OpenGL state
    mover_t mover_GLstate;
    init_mover(&mover_GLstate, this->screen_width(), this->screen_height());

    // Initialize font OpenGL state
    font_t font_state;
    init_font(&font_state, this->screen_width(), this->screen_height());

    // Init container OpenGL state
    container_t container_state;
    init_container(&container_state, this->screen_width(), this->screen_height());

    // Initialize RGB Light if present
    #if defined BLINK1
    rgb_light_t light_state;
    init_rgb_light(&light_state, 255, 0, 0);
    #endif

    // Broadcast pixels ratio
    short pixel_dims[2];
    pixel_dims[0] = (short)this->screen_width();
    pixel_dims[1] = (short)this->screen_height();
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
 
    // Recv simulation world dimensions from global rank 1
    float sim_dims[3];
    MPI_Recv(sim_dims, 3, MPI_FLOAT, 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    this->sim_width = sim_dims[0];
    this->sim_depth = sim_dims[1];
    this->sim_height = sim_dims[2];
    // Receive number of global particles
    int max_particles;
    MPI_Recv(&max_particles, 1, MPI_INT, 1, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Gatherv initial tunable parameters values
    int *param_counts = (int*)malloc(this->num_compute_procs+1 * sizeof(int));
    int *param_displs = (int*)malloc(this->num_compute_procs+1 * sizeof(int));
    for(i=0; i<this->num_compute_procs+1; i++) {
        param_counts[i] = i?1:0; // will not receive from rank 0
        param_displs[i] = i?i-1:0; // rank i will reside in params[i-1]
    }
    // Initial gather from compute nodes to renderer
    MPI_Gatherv(MPI_IN_PLACE, 0, TunableParamtype, this->tunable_param_structs.data(), param_counts, param_displs, TunableParamtype, 0, MPI_COMM_WORLD);
    this->param_struct_to_class();

    // Allocate particle receive array
    int num_coords = 3;
    short *particle_coords = (short*)malloc(num_coords * max_particles*sizeof(short));

    // Allocate points array(position + color)
    int point_size = 6 * sizeof(float);
    float *points = (float*)malloc(point_size*max_particles);

    // Allocate mover point array(position + color)
    float mover_center[3];
    float mover_color[3];

    // Number of coordinates received from each proc
    int *particle_coordinate_counts = (int*)malloc(num_compute_procs * sizeof(int));
    // Keep track of order in which particles received
    int *particle_coordinate_ranks = (int*)malloc(num_compute_procs * sizeof(int));

    // Create color index, equally spaced around HSV
    float *colors_by_rank = (float*)malloc(3*this->num_compute_procs*sizeof(float));
    float angle_space = 0.5f/(float)this->num_compute_procs;
    float HSV[3];
    for(i=0; i<this->num_compute_procs; i++)
    {
        if(i%2)
            HSV[0] = angle_space*i;
        else
            HSV[0] = angle_space*i + 0.5f;
        HSV[1] = 1.0f;
        HSV[2] = 0.8f;
        hsv_to_rgb(HSV, colors_by_rank+3*i);
    }
 
    #if defined BLINK1
    MPI_Bcast(colors_by_rank, 3*this->num_compute_procs, MPI_FLOAT, 0, MPI_COMM_WORLD);
    #endif

    int num_coords_rank;
    int current_rank, num_parts;

    int frames_per_fps = 30;
    int frames_per_check = 1;
    int num_steps = 0;
    double current_time;
    double wall_time = MPI_Wtime();
    float fps=0.0f;

    // Setup MPI requests used to gather particle coordinates
    MPI_Request coord_reqs[num_compute_procs];
    int src, coords_recvd;
    float gl_x, gl_y;
    // Particle radius in pixels
    float particle_diameter_pixels = gl_state.screen_width * 0.0425;
    float liquid_particle_diameter_pixels = gl_state.screen_width * 0.015;

    MPI_Status status;

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
            this->tunable_parameters->kill_sim();
            this->param_class_to_struct();
            // Send kill paramaters to compute nodes
            MPI_Scatterv(node_params, param_counts, param_displs, TunableParamtype, MPI_IN_PLACE, 0, TunableParamtype, 0, MPI_COMM_WORLD);
            break;
        }    

        // Check for user keyboard/mouse input
        if(this->pause) {
            while(this->pause)
                this->glfw->check_user_input();
        }
        else
            this->glfw->check_user_input();

        // Check if inactive
//        if(!this->input_is_active())
//            this->update_inactive_state();

        // Update struct params with class values
        this->param_class_to_struct();

        // Send updated paramaters to compute nodes
        MPI_Scatterv(node_params, param_counts, param_displs, TunableParamtype, MPI_IN_PLACE, 0, TunableParamtype, 0, MPI_COMM_WORLD);

            // Retrieve all particle coordinates (x,y,z)
  	    // Potentially probe is expensive? Could just allocated num_compute_procs*num_particles_global and async recv
	    // OR do synchronous recv...very likely that synchronous receive is as fast as anything else
	    coords_recvd = 0;
	    for(i=0; i<this->num_compute_procs; i++) {
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
            this->tunable_parameters->check_partition_left(particle_coordinate_counts, coords_recvd);

        // Clear background
        glClearColor(0.15, 0.15, 0.15, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw container image
        render_container(&container_state);

        // update mover
        this->sim_to_opengl(render_state.master_params[0].mover_center_x,
                                     render_state.master_params[0].mover_center_y,
                                     render_state.master_params[0].mover_center_z,
                                     &mover_center[0], &mover_center[1], &mover_center[2]);

        float mover_radius = this->mover_radius/this->sim_width * 1.0f;
        mover_color[0] = 1.0f;
        mover_color[1] = 0.0f;
        mover_color[2] = 0.0f;

        render_all_text(&font_state, &render_state, fps);

        // Wait for all coordinates to be received
        MPI_Waitall(num_compute_procs, coord_reqs, MPI_STATUSES_IGNORE);

        // Create points array (x,y,r,g,b)
        i = 0;
        current_rank = particle_coordinate_ranks[i];
        // j == coordinate pair
        for(j=0, num_parts=1; j<coords_recvd/3; j++, num_parts++) {
             // Check if we are processing a new rank's particles
             if ( num_parts > particle_coordinate_counts[current_rank]/3){
                current_rank =  particle_coordinate_ranks[++i];
                num_parts = 1;
                // Find next rank with particles if current_rank has 0 particles
                while(!particle_coordinate_counts[current_rank])
                    current_rank = particle_coordinate_ranks[++i];
            }
            points[j*6]   = particle_coords[j*3]/(float)SHRT_MAX;
            points[j*6+1] = particle_coords[j*3+1]/(float)SHRT_MAX;
            points[j*6+2] = particle_coords[j*3+2]/(float)SHRT_MAX;
            points[j*6+3] = colors_by_rank[3*current_rank];
            points[j*6+4] = colors_by_rank[3*current_rank+1];
            points[j*6+5] = colors_by_rank[3*current_rank+2];
        }

        render_particles(points, particle_diameter_pixels, coords_recvd/3, &particle_GLstate);

        // Render over particles to hide penetration
        render_mover(mover_center, mover_radius, mover_color, &mover_GLstate);

        // Swap front/back buffers
        this->glfw->swap_ogl(&gl_state);

        update_view(&world_GLstate);

        num_steps++;
    }

    #if defined BLINK1
    shutdown_rgb_light(&light_state);
    #endif

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
void Renderer::opengl_to_sim(float x, float y, float z, float *sim_x, float *sim_y, float *sim_z)
{
    float half_scale = this->sim_width*0.5f;

    *sim_x = x*half_scale + half_scale;
    *sim_y = y*half_scale + half_scale;
    *sim_z = z*half_scale + half_scale;
}

// Translate between simulation coordinates, origin bottom left, and opengl -1,1 center of screen coordinates
void Renderer::sim_to_opengl(render_t *render_state, float x, float y, float z, float *gl_x, float *gl_y, float *gl_z)
{
    float half_scale = this->sim_width*0.5f;

    *gl_x = x/half_scale - 1.0f;
    *gl_y = y/half_scale - 1.0f;
    *gl_z = z/half_scale - 1.0f;
}

// Set time of last user input
void Renderer::set_activity_time()
{
    this->last_activity_time = MPI_Wtime();
}

// Return if simulation has active input or not
bool Renderer::input_is_active()
{
    double time_since_active = MPI_Wtime() - this->last_activity_time;
    return time_since_active < 120;
}

void Renderer::enable_view_controls()
{
    this->view_controls = true;
}

void Renderer::disable_view_controls()
{
    this->view_controls = false;
}

void Renderer::toggle_pause()
{
    this->pause = !state->pause;
}

void Renderer::set_view_angle(float x_pos, float y_pos)
{
    world_t *world_state = this->world;
    const float max_degrees_rotate = 40.0f;

    // x_pos,y_pos is [-1, 1]
    // Angle is [-max_degrees, max_degrees]
    // This is the angle relative to the intial orientation at the start of view mode
    float degrees_yaw = max_degrees*x_pos;
    float degrees_pitch = max_degrees*y_pos;
    rotate_camera_yaw_pitch(world_state, degrees_yaw, degrees_pitch);

    // Update view
    update_view(world_state);
}

void Renderer::move_in_view()
{
    float dx = 0.15;
    world_t *world_state = this->world;
    move_along_view(world_state, dx);
    update_view(world_state);
}

void Renderer::move_out_view()
{
    float dx = -0.15;
    world_t *world_state = state->world;
    move_along_view(world_state, dx);
    update_view(world_state);
}

void Renderer::zoom_in_view()
{
    float dzoom = 0.07;
    world_t *world_state = state->world;
    zoom_view(world_state, dzoom);
    update_view(world_state);
}

void Renderer::zoom_out_view()
{
    float dzoom = -0.07;
    world_t *world_state = state->world;
    zoom_view(world_state, dzoom);
    update_view(world_state);
}
