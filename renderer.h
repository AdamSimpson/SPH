#ifndef fluid_renderer_h
#define fluid_renderer_h

#include "fluid.h"

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

typedef struct {
    parameters selected_parameter;
    tunable_parameters *node_params; // Holds all nodes paramters including start/end lengths
    tunable_parameters master_params; // Holds parameters shared by all nodes
    int num_compute_procs;
} RENDER_T;

void start_renderer();
void increase_parameter(RENDER_T *render_state);
void decrease_parameter(RENDER_T *render_state);
void move_parameter_up(RENDER_T *render_state);
void move_parameter_down(RENDER_T *render_state);
void pixel_to_sim(float *world_dims, float x, float y, float *sim_x, float *sim_y);
void sim_to_opengl(float *world_dims, float x, float y, float *gl_x, float *gl_y);
void sim_pixels_to_opengl(short *pixel_dims, short x, short y, float *gl_x, float *gl_y);
void decrease_gravity(RENDER_T *render_state);
void increase_gravity(RENDER_T *render_state);

#endif
