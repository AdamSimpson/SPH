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

void start_renderer()
{
    // Setup initial OpenGL state
    gl_t gl_state;
    memset(&gl_state, 0, sizeof(gl_t));

    // Start OpenGL
    render_t render_state;
    init_ogl(&gl_state, &render_state);

    // Initialize particles OpenGL state
    particles_t particle_GLstate;
    init_particles(&particle_GLstate, gl_state.screen_width, gl_state.screen_height);

    // Initialize mover OpenGL state
    mover_t mover_GLstate;
    init_mover(&mover_GLstate);

    // Initialize font atlas
    font_t font_state;
    init_font(&font_state, gl_state.screen_width, gl_state.screen_height);

    // Initialize background OpenGL state
    background_t background_state;
    init_background(&background_state);

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
 
    // Recv world dimensions from global rank 1
    float world_dims[2];
    MPI_Recv(world_dims, 2, MPI_FLOAT, 1, 8, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Receive number of global particles
    int max_particles;
    MPI_Recv(&max_particles, 1, MPI_INT, 1, 9, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Calculate world unit to pixel
    float world_to_pix_scale = gl_state.screen_width/world_dims[0];

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
    float *points = (float*)malloc(point_size*max_particles);

    // Allocate mover point array(position + color)
    float mover_center[2];
    float mover_color[3];

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
 
    int num_coords_rank;
    int current_rank, num_parts;
    float mouse_x, mouse_y, mouse_x_scaled, mouse_y_scaled;
    float mover_gl_dims[2];

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
    float particle_diameter_pixels = gl_state.screen_width * 0.0125;

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
            for(i=0; i<render_state.num_compute_procs; i++)
                render_state.node_params[i].kill_sim = true;
            // Send kill paramaters to compute nodes
            MPI_Scatterv(node_params, param_counts, param_displs, TunableParamtype, MPI_IN_PLACE, 0, TunableParamtype, 0, MPI_COMM_WORLD);
            break;
        }    

        // Update mover position
        get_mouse(&mouse_x, &mouse_y, &gl_state);
        pixel_to_sim(world_dims, mouse_x, mouse_y, &mouse_x_scaled, &mouse_y_scaled);

        // Mover position in sim coordinates
        for(i=0; i<render_state.num_compute_procs; i++) {
            render_state.master_params[i].mover_center_x = mouse_x_scaled;
            render_state.master_params[i].mover_center_y = mouse_y_scaled;
        }

        // Get keyboard key press
        // process appropriately
        check_key_press(&gl_state);

        // Update node params with master param values
        update_node_params(&render_state);

        // Send updated paramaters to compute nodes
        MPI_Scatterv(node_params, param_counts, param_displs, TunableParamtype, MPI_IN_PLACE, 0, TunableParamtype, 0, MPI_COMM_WORLD);

        // Retrieve all particle coordinates (x,y)
  	    // Potentially probe is expensive? Could just allocated num_compute_procs*num_particles_global and async recv
	    // OR do synchronous recv...very likely that synchronous receive is as fast as anything else
	    coords_recvd = 0;
	    for(i=0; i<num_compute_procs; i++) {
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
        check_partition_right(&render_state, particle_coordinate_counts, coords_recvd);

        // Clear background
        glClear(GL_COLOR_BUFFER_BIT);

        // Draw background image
        draw_background(&background_state);

        // update mover
        sim_to_opengl(world_dims, mouse_x_scaled, mouse_y_scaled, &gl_x, &gl_y);
        mover_center[0] = gl_x;
        mover_center[1] = gl_y;
        mover_color[0] = 1.0f;
        mover_color[1] = 0.0f;
        mover_color[2] = 0.0f;

        // Mover bounding rectangle half width/height lengths in ogl system
        // Subtract off particle diamter so no particle/mover penetration
        mover_gl_dims[0] = render_state.master_params[0].mover_width/(world_dims[0]*0.5f) - particle_diameter_pixels/(gl_state.screen_width*0.5f) ;
        mover_gl_dims[1] = render_state.master_params[0].mover_height/(world_dims[1]*0.5f) - particle_diameter_pixels/(gl_state.screen_height*0.5f);

        render_all_text(&font_state, &render_state, fps);
        //printf("fps: %f\n", fps);

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

    exit_ogl(&gl_state);
}

// Move selected parameter up
void move_parameter_up(render_t *render_state)
{
    if(render_state->selected_parameter == MIN)
        render_state->selected_parameter = MAX;
    else
	render_state->selected_parameter--;
}

// Move selected parameter down
void move_parameter_down(render_t *render_state) 
{
    if(render_state->selected_parameter == MAX)
        render_state->selected_parameter = MIN;
    else
        render_state->selected_parameter++;
}

void increase_parameter(render_t *render_state)
{
    switch(render_state->selected_parameter) {
        case GRAVITY:
	        increase_gravity(render_state);
	        break;
        case VISCOSITY:
            increase_viscosity(render_state);
            break;
        case DENSITY:
            increase_density(render_state);
            break;
        case PRESSURE:
            increase_pressure(render_state);
            break;
         case ELASTICITY:
            increase_elasticity(render_state);
            break;
    }
}

void decrease_parameter(render_t *render_state)
{
    switch(render_state->selected_parameter) {
        case GRAVITY:
            decrease_gravity(render_state);
            break;
        case VISCOSITY:
            decrease_viscosity(render_state);
            break;
        case DENSITY:
            decrease_density(render_state);
            break;
        case PRESSURE:
            decrease_pressure(render_state);
            break;
        case ELASTICITY:
            decrease_elasticity(render_state);
            break;
    }
}

// Increase gravity parameter
void increase_gravity(render_t *render_state)
{
    static const float max_grav = -9.0f;
    if(render_state->master_params[0].g <= max_grav)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->master_params[i].g -= 1.0f;
}

// Decreate gravity parameter
void decrease_gravity(render_t *render_state)
{
    static const float min_grav = 9.0f;
    if(render_state->master_params[0].g >= min_grav)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->master_params[i].g += 1.0f;
}

// Increase density parameter
void increase_density(render_t *render_state)
{
    static const float max_dens = 150.0f;
    if(render_state->master_params[0].rest_density >= max_dens)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->master_params[i].rest_density += 5.0f;
}

// Decreate gravity parameter
void decrease_density(render_t *render_state)
{
    static const float min_dens = 0.0f;
    if(render_state->master_params[0].rest_density <= min_dens)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->master_params[i].rest_density -= 5.0f;
}

// Increase viscosity parameter
void increase_viscosity(render_t *render_state)
{
    static const float max_viscosity = 100.0f;
    float viscosity = render_state->master_params[0].sigma;

    if(viscosity > max_viscosity)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].sigma += 5.0f;
        render_state->master_params[i].beta += 0.5f;
    }
}

// Decreate viscosity parameter
void decrease_viscosity(render_t *render_state)
{
    static const float min_viscosity = 0.0f;
    float viscosity = render_state->master_params[0].sigma;

    if(viscosity <= min_viscosity)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].sigma -= 5.0f;
        render_state->master_params[i].beta -= 0.5f;
    }
}

// Increase pressure parameter
void increase_pressure(render_t *render_state)
{
    static const float max_pressure = 1.0f;
    float pressure = render_state->master_params[0].k;
    float pressure_near = render_state->master_params[0].k_near;

    if(pressure > max_pressure)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].k += 0.1f;
        render_state->master_params[i].k_near += 0.5f;
    }
}

// Decreate pressure parameter
void decrease_pressure(render_t *render_state)
{
    static const float min_pressure = 0.0f;
    float pressure = render_state->master_params[0].k;
    float pressure_near = render_state->master_params[0].k_near;

    if(pressure <= min_pressure)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].k -= 0.1f;
        render_state->master_params[i].k_near -= 0.5f;
    }
}

// Increase elasticity parameter
void increase_elasticity(render_t *render_state)
{
    static const float max_elast = 200.0f;
    if(render_state->master_params[0].k_spring > max_elast)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->master_params[i].k_spring += 5.0f;
}

// Decreate elasticity parameter
void decrease_elasticity(render_t *render_state)
{
    static const float min_elast = -50.0f;
    if(render_state->master_params[0].k_spring < min_elast)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++)
        render_state->master_params[i].k_spring -= 5.0f;
}

// Increase mover x dimension
void increase_mover_width(render_t *render_state)
{
    // Maximum width of mover
    static const float max_width = 4.0f;

    if(render_state->master_params[0].mover_width > max_width)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        // Increase sphere
        if(render_state->master_params[i].mover_type == SPHERE_MOVER) {
            render_state->master_params[i].mover_height += 0.2f;
            render_state->master_params[i].mover_width += 0.2f;
        }
        else if(render_state->master_params[i].mover_type == RECTANGLE_MOVER) {
            render_state->master_params[i].mover_width += 0.2f;
        }
    }

}

// Decrease mover x dimension
void decrease_mover_width(render_t *render_state)
{
    // Minimum width of mover
    static const float min_width = 1.0f;

    if(render_state->master_params[0].mover_width - min_width < 0.001f)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        // Decrease sphere radius
        if(render_state->master_params[i].mover_type == SPHERE_MOVER) {
            render_state->master_params[i].mover_height -= 0.2f;
            render_state->master_params[i].mover_width -= 0.2f;
        }
        // Decrease rectangle width
        else if(render_state->master_params[i].mover_type == RECTANGLE_MOVER) {
            render_state->master_params[i].mover_width -= 0.2f;
        }
    }
}

// Increase mover y dimension
void increase_mover_height(render_t *render_state)
{
    // Maximum height of mover
    static const float max_height = 4.0f;

    if(render_state->master_params[0].mover_height > max_height)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        // Increase sphere radius
        if(render_state->master_params[i].mover_type == SPHERE_MOVER) {
            render_state->master_params[i].mover_height += 0.2f;
            render_state->master_params[i].mover_width += 0.2f;
        }
        // Increase rectangle height
        else if(render_state->master_params[i].mover_type == RECTANGLE_MOVER) {
            render_state->master_params[i].mover_height += 0.2f;
        }
    }

}

// Decrease mover y dimension
void decrease_mover_height(render_t *render_state)
{
    // Minimum height of mover
    static const float min_height = 1.0f;

    if(render_state->master_params[0].mover_height - min_height < 0.001f)
        return;

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        // Decrease sphere radius
        if(render_state->master_params[i].mover_type == SPHERE_MOVER) {
            render_state->master_params[i].mover_height -= 0.2f;
            render_state->master_params[i].mover_width -= 0.2f;
        }
        // Decrease rectangle height
        else if(render_state->master_params[i].mover_type == RECTANGLE_MOVER) {
            render_state->master_params[i].mover_height -= 0.2f;
        }
    }
}

// Preset fluids based upon xbox controller buttons

// Standard fluid
void set_fluid_x(render_t *render_state)
{
    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].g = 6.0f;
        render_state->master_params[i].k = 0.2f;
        render_state->master_params[i].k_near = 6.0f;
        render_state->master_params[i].k_spring = 10.0f;
        render_state->master_params[i].sigma = 5.0f;
        render_state->master_params[i].beta = 0.5f;
        render_state->master_params[i].rest_density = 30.0f;  
    }
}

// super goo
void set_fluid_y(render_t *render_state)
{
    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].g = 6.0f;
        render_state->master_params[i].k = 0.1f;
        render_state->master_params[i].k_near = 3.0f;
        render_state->master_params[i].k_spring = -30.0f;
        render_state->master_params[i].sigma = 100.0f;
        render_state->master_params[i].beta = 10.0f;
        render_state->master_params[i].rest_density = 30.0f;
    }
}

// zero g high surface tension
void set_fluid_a(render_t *render_state)
{
    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].g = 0.0f;
        render_state->master_params[i].k = 0.2f;
        render_state->master_params[i].k_near = 6.0f;
        render_state->master_params[i].k_spring = 10.0f;
        render_state->master_params[i].sigma = 20.0f;
        render_state->master_params[i].beta = 2.0f;
        render_state->master_params[i].rest_density = 55.0f;
    }
}

// spring gas fluid...thing
void set_fluid_b(render_t *render_state)
{
    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].g = 6.0f;
        render_state->master_params[i].k = 0.0f;
        render_state->master_params[i].k_near = 0.0f;
        render_state->master_params[i].k_spring = 115.0f;
        render_state->master_params[i].sigma = 20.0f;
        render_state->master_params[i].beta = 2.0f;
        render_state->master_params[i].rest_density = 0.0f;
    }
}

// Translate between pixel coordinates with origin at screen center
// to simulation coordinates
void pixel_to_sim(float *world_dims, float x, float y, float *sim_x, float *sim_y)
{
    float half_width = world_dims[0]*0.5f;
    float half_height = world_dims[1]*0.5f;

    *sim_x = x*half_width + half_width;
    *sim_y = y*half_height + half_height;
}

// Translate between simulation coordinates, origin bottom left, and opengl -1,1 center of screen coordinates
void sim_to_opengl(float *world_dims, float x, float y, float *gl_x, float *gl_y)
{
    float half_width = world_dims[0]*0.5f;
    float half_height = world_dims[1]*0.5f;

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
// Check from left to right
void check_partition_right(render_t *render_state, int *particle_counts, int total_particles)
{
    int rank, diff;
    float h, dx, length, length_right;

    // Particles per proc if evenly divided
    int even_particles = total_particles/render_state->num_compute_procs_active;
    int max_diff = even_particles/15.0f;

    // Fixed distance to move partition is 1/2 smoothing radius
    h = render_state->master_params[0].smoothing_radius;
    dx = h*0.5;

    tunable_parameters *master_params = render_state->master_params;

    for(rank=0; rank<(render_state->num_compute_procs_active-1); rank++)
    {
        length =  master_params[rank].node_end_x - master_params[rank].node_start_x; 
        length_right =  master_params[rank+1].node_end_x - master_params[rank+1].node_start_x; 
        diff = particle_counts[rank] - even_particles;// particle_counts[rank+1];

        // current rank has too many particles
        if( diff > max_diff && length > 4*h) {
            master_params[rank].node_end_x -= dx;
            master_params[rank+1].node_start_x = master_params[rank].node_end_x;
	    rank++;
        }
        // current rank has too few particles
        else if (diff < -max_diff && length_right > 4*h) {
            master_params[rank].node_end_x += dx;
            master_params[rank+1].node_start_x = master_params[rank].node_end_x;
	    rank++;
        }    
    }
}

// Checks for a balanced number of particles on each compute node
// If unbalanced the partition between nodes will change 
// Check from right to left
void check_partition_left(render_t *render_state, int *particle_counts, int total_particles)
{
    int rank, diff;
    float h, dx, length, length_left;   

    // Particles per proc if evenly divided
    int even_particles = total_particles/render_state->num_compute_procs_active;
    int max_diff = even_particles/10.0f;

    // Fixed distance to move partition is 1/2 smoothing radius
    h = render_state->master_params[0].smoothing_radius;
    dx = h*0.125;

    tunable_parameters *master_params = render_state->master_params;

    for(rank=render_state->num_compute_procs_active; rank-- > 1; )
    {
        length =  master_params[rank].node_end_x - master_params[rank].node_start_x;
        length_left =  master_params[rank-1].node_end_x - master_params[rank-1].node_start_x;
        diff = particle_counts[rank] - even_particles;//particle_counts[rank-1];

        // current rank has too many particles
        if( diff > max_diff && length > 4*h) {
            master_params[rank].node_start_x += dx;
            master_params[rank-1].node_end_x = master_params[rank].node_start_x;
        }
        // current rank has too few particles
        else if (diff < -max_diff && length_left > 4*h) {
            master_params[rank].node_start_x -= dx;
            master_params[rank-1].node_end_x = master_params[rank].node_start_x;
        }

    }
}

// Set last partition to be outside of simulation bounds
// Effectively removing it from the simulation
void remove_partition(render_t *render_state)
{
    if(render_state->num_compute_procs_active == 1) 
	return;

    int num_compute_procs_active = render_state->num_compute_procs_active;

    // Set new end position of last active proc to end of simulation
    render_state->master_params[num_compute_procs_active-2].node_end_x = render_state->master_params[num_compute_procs_active-1].node_end_x;

    // Send start and end x out of sim bounds
    float position = render_state->master_params[num_compute_procs_active-1].node_end_x + 1.0; // +1.0 ensures it's out of the simulation bounds
    render_state->master_params[num_compute_procs_active-1].node_start_x = position;
    render_state->master_params[num_compute_procs_active-1].node_end_x = position;

    render_state->num_compute_procs_active -= 1;
}

// Add on partition to right side that has been removed
void add_partition(render_t *render_state)
{
    if(render_state->num_compute_procs_active == render_state->num_compute_procs)
	return;

    // Length of currently last partiion
    int num_compute_procs_active = render_state->num_compute_procs_active;
    float length = render_state->master_params[num_compute_procs_active-1].node_end_x - render_state->master_params[num_compute_procs_active-1].node_start_x;
    float h = render_state->master_params[0].smoothing_radius;

    // If the last partition is too small we can't split it and another
    if(length < 6*h)
	return;

    // Set end of added partition to current end location
    render_state->master_params[num_compute_procs_active].node_end_x = render_state->master_params[num_compute_procs_active-1].node_end_x;
    
    // Divide the current last partition in half
    float new_x = render_state->master_params[num_compute_procs_active-1].node_start_x + length*0.5;
    render_state->master_params[num_compute_procs_active-1].node_end_x = new_x;
    render_state->master_params[num_compute_procs_active].node_start_x = new_x;

    render_state->num_compute_procs_active += 1;
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
