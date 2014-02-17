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
    tunable_parameters *master_params; // Holds parameters shared by all nodes
    int num_compute_procs;
    int num_compute_procs_active; // Number of nodes participating in simulation, user may "remove" nodes at runtime
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
void decrease_viscosity(RENDER_T *render_state);
void increase_viscosity(RENDER_T *render_state);
void decrease_density(RENDER_T *render_state);
void increase_density(RENDER_T *render_state);
void decrease_pressure(RENDER_T *render_state);
void increase_pressure(RENDER_T *render_state);
void decrease_elasticity(RENDER_T *render_state);
void increase_elasticity(RENDER_T *render_state);
void update_node_params(RENDER_T *render_state);
void checkPartitions(RENDER_T *render_state, int *particle_counts, int total_particles);
void add_partition(RENDER_T *render_state);
void remove_partition(RENDER_T *render_state);
void hsv_to_rgb(float* hsv, float *rgb);
void check_partition_right(RENDER_T *render_state, int *particle_counts, int total_particles);
void check_partition_left(RENDER_T *render_state, int *particle_counts, int total_particles);

#endif
