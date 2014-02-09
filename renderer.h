#ifndef fluid_renderer_h
#define fluid_renderer_h

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

void start_renderer();
void increase_parameter();
void decrease_parameter();
void move_parameter_up();
void move_parameter_down();
void pixel_to_sim(float *world_dims, float x, float y, float *sim_x, float *sim_y);
void sim_to_opengl(float *world_dims, float x, float y, float *gl_x, float *gl_y);

#endif
