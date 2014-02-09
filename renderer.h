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
    param *params;
    int num_params;

    // "Master" parameters
    float rest_density;
    float spacing_particle;
    float smoothing_radius;
    float g;
    float k; // Pressure constant
    float k_near; // Near pressure constant
    float k_spring; // Spring constant
    float sigma; // linear velocity viscocity term
    float beta;  // quadratic velocity viscocity term
    float time_step;
    float node_start_x; // left x position of node partition
    float node_end_x;   // right x position of node partition
    float mover_center_x;
    float mover_center_y;
    float mover_radius;
} RENDER_T;

void start_renderer();
void increase_parameter(RENDER_T *render_state);
void decrease_parameter(RENDER_T *render_state);
void move_parameter_up(RENDER_T *render_state);
void move_parameter_down(RENDER_T *render_state);
void pixel_to_sim(float *world_dims, float x, float y, float *sim_x, float *sim_y);
void sim_to_opengl(float *world_dims, float x, float y, float *gl_x, float *gl_y);
void decrease_gravity(RENDER_T *render_state);
void increase_gravity(RENDER_T *render_state);

#endif
