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

typedef struct render_t {
    parameters selected_parameter;
    tunable_parameters *node_params; // Holds all nodes paramters including start/end lengths
    tunable_parameters *master_params; // Holds parameters shared by all nodes
    int num_compute_procs;
    int num_compute_procs_active; // Number of nodes participating in simulation, user may "remove" nodes at runtime
} render_t;

void start_renderer();
void increase_parameter(render_t *render_state);
void decrease_parameter(render_t *render_state);
void move_parameter_up(render_t *render_state);
void move_parameter_down(render_t *render_state);
void pixel_to_sim(float *world_dims, float x, float y, float *sim_x, float *sim_y);
void sim_to_opengl(float *world_dims, float x, float y, float *gl_x, float *gl_y);
void sim_pixels_to_opengl(short *pixel_dims, short x, short y, float *gl_x, float *gl_y);
void decrease_gravity(render_t *render_state);
void increase_gravity(render_t *render_state);
void decrease_viscosity(render_t *render_state);
void increase_viscosity(render_t *render_state);
void decrease_density(render_t *render_state);
void increase_density(render_t *render_state);
void decrease_pressure(render_t *render_state);
void increase_pressure(render_t *render_state);
void decrease_elasticity(render_t *render_state);
void increase_elasticity(render_t *render_state);
void update_node_params(render_t *render_state);
void checkPartitions(render_t *render_state, int *particle_counts, int total_particles);
void add_partition(render_t *render_state);
void remove_partition(render_t *render_state);
void hsv_to_rgb(float* hsv, float *rgb);
void check_partition_right(render_t *render_state, int *particle_counts, int total_particles);
void check_partition_left(render_t *render_state, int *particle_counts, int total_particles);
void decrease_mover_height(render_t *render_state);
void increase_mover_height(render_t *render_state);
void decrease_mover_width(render_t *render_state);
void increase_mover_width(render_t *render_state);
void set_fluid_x(render_t *render_state);
void set_fluid_y(render_t *render_state);
void set_fluid_a(render_t *render_state);
void set_fluid_b(render_t *render_state);

#endif
