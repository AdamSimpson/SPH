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

#include "controls.h"
#include "renderer.h"

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

// Set center of mover, input is openGL coordinates
void set_mover_gl_center(render_t *render_state, float ogl_x, float ogl_y)
{
    float sim_x, sim_y;
    opengl_to_sim(render_state, ogl_x, ogl_y, &sim_x, &sim_y);

    int i;
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].mover_center_x = sim_x;
        render_state->master_params[i].mover_center_z = sim_y;
        render_state->master_params[i].mover_center_y = 5.0f;
    }  
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
        render_state->master_params[i].mover_height += 0.2f;
        render_state->master_params[i].mover_width += 0.2f;
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
        render_state->master_params[i].mover_height -= 0.2f;
        render_state->master_params[i].mover_width -= 0.2f;
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
        render_state->master_params[i].mover_height += 0.2f;
        render_state->master_params[i].mover_width += 0.2f;
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
        render_state->master_params[i].mover_height -= 0.2f;
        render_state->master_params[i].mover_width -= 0.2f;
    }
}

// Reset the mover width and height
void reset_mover_size(render_t *render_state) {
    int i;     
    for(i=0; i<render_state->num_compute_procs; i++) {
        render_state->master_params[i].mover_height = 2.0f;
        render_state->master_params[i].mover_width = 2.0f;
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

// Set last partition to be outside of simulation bounds
// Effectively removing it from the simulation
void remove_partition(render_t *render_state)
{
    if(render_state->num_compute_procs_active == 1) 
	return;

    int num_compute_procs_active = render_state->num_compute_procs_active;

    int removed_rank = num_compute_procs_active-1;

    // Set new end position of last active proc to end of simulation
    render_state->master_params[removed_rank-1].node_end_x = render_state->master_params[removed_rank].node_end_x;

    // Send start and end x out of sim bounds
    float position = render_state->master_params[removed_rank].node_end_x + 1.0; // +1.0 ensures it's out of the simulation bounds
    render_state->master_params[removed_rank].node_start_x = position;
    render_state->master_params[removed_rank].node_end_x = position;

    // Set active to false for removed rank
    render_state->master_params[removed_rank].active = false;

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
    if(length < 2.5*h)
	return;

    // Set end of added partition to current end location
    render_state->master_params[num_compute_procs_active].node_end_x = render_state->master_params[num_compute_procs_active-1].node_end_x;
    
    // Divide the current last partition in half
    float new_x = render_state->master_params[num_compute_procs_active-1].node_start_x + length*0.5;
    render_state->master_params[num_compute_procs_active-1].node_end_x = new_x;
    render_state->master_params[num_compute_procs_active].node_start_x = new_x;

    // Set active to true for added rank
    render_state->master_params[num_compute_procs_active].active = true;

    render_state->num_compute_procs_active += 1;
}

void toggle_dividers(render_t *state)
{
    state->show_dividers = !state->show_dividers;
}

void toggle_pause(render_t *state)
{
    state->pause = !state->pause;
}

void toggle_liquid(render_t *state)
{
    state->liquid = !state->liquid;
}
