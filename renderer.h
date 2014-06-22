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

#ifndef fluid_renderer_h
#define fluid_renderer_h

#include "fluid.h"
#include "stdbool.h"

// enum of displayed parameter values
typedef enum {
    MIN = 0,
    GRAVITY = MIN,
    VISCOSITY,
    DENSITY,
    PRESSURE,
    ELASTICITY,
    MAX = ELASTICITY
 } parameters;

typedef struct render_t {
    float sim_width;
    float sim_height;
    float screen_width;
    float screen_height;
    parameters selected_parameter;
    tunable_parameters *node_params; // Holds all nodes paramters including start/end lengths
    tunable_parameters *master_params; // Holds parameters shared by all nodes
    int num_compute_procs;
    int num_compute_procs_active; // Number of nodes participating in simulation, user may "remove" nodes at runtime
    bool show_dividers;
    bool pause;
    double last_activity_time; // Used to determine if simulation is being used or not
} render_t;

void start_renderer();
void opengl_to_sim(render_t *render_state, float x, float y, float *sim_x, float *sim_y);
void sim_to_opengl(render_t *render_state, float x, float y, float *gl_x, float *gl_y);
void update_node_params(render_t *render_state);
void checkPartitions(render_t *render_state, int *particle_counts, int total_particles);
void hsv_to_rgb(float* hsv, float *rgb);
void check_partition_left(render_t *render_state, int *particle_counts, int total_particles);
void set_activity_time(render_t *render_state);
bool input_is_active(render_t *render_state);
void update_inactive_state(render_t *render_state);

#endif
