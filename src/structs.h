#ifndef __STRUCTS_H
#define __STRUCTS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "stdbool.h"

typedef unsigned int uint;

typedef struct gl_t {
    int screen_width;
    int screen_height;

    GLFWwindow* window;
} gl_t;

// Standard fluid particle paramaters
typedef struct FLUID_PARTICLES_T {
    float *x_star;
    float *y_star;
    float *z_star;
    float *x;
    float *y;
    float *z;
    float *v_x;
    float *v_y;
    float *v_z;
    float *dp_x;
    float *dp_y;
    float *dp_z;
    float *density;
    float *lambda;
    uint *id; // Id is 'local' index within the fluid particle pointer array
} fluid_particles_t;

// Particles that are within 2*h distance of node edge
typedef struct EDGE_T {
    int max_edge_particles;
    uint *edge_indices_left;
    uint *edge_indices_right;
    int number_edge_particles_left;
    int number_edge_particles_right;
} edge_t;

// Particles that have left the node
typedef struct OOB_T {
    int max_oob_particles;
    uint *oob_index_indices_left; // Indicies in particle index array for particles traveling left
    uint *oob_index_indices_right;
    int number_oob_particles_left; // Number of OOB sending left
    int number_oob_particles_right;
    uint *vacant_indices; // Indicies in global particle array that are vacant
    int number_vacancies;
} oob_t;

typedef struct AABB_T {
    float min_x;
    float max_x;
    float min_y;
    float max_y;
    float min_z;
    float max_z;
} AABB_t; //Axis aligned bounding box

// Bucket to hold each particles nearest neighbors
typedef struct NEIGHBOR_T{
    uint *fluid_neighbors; // Index in global particle array of neighbor
    int number_fluid_neighbors;
} neighbor_t;

typedef struct NEIGHBOR_GRID_T {
    float spacing;  // Spacing between buckets
    uint size_x; // Number of buckets in x
    uint size_y; // Number of buckets in y
    uint size_z;
    uint *start_indices; // Start index for hash values
    uint *end_indices;   // End index for hash values
    uint *hash_values; // Array of hash values
    uint *particle_ids; // Array of particle id's
    uint max_neighbors; // Maximum neighbors allowed for each particle
    neighbor_t *neighbors; // Particle neighbor buckets
} neighbor_grid_t;

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

// These parameters are tunable by the render node
typedef struct TUNABLE_PARAMETERS_T {
    float rest_density;
    float smoothing_radius;
    float g;
    float k;
    float k_near;
    float k_spring;
    float sigma;
    float beta;
    float time_step;
    float node_start_x;
    float node_end_x;
    float mover_center_x;
    float mover_center_y;
    float mover_center_z;
    float mover_width;
    float mover_height;
    char mover_type;
    char kill_sim;
    char active;
} tunable_parameters_t;

// Full parameters struct for simulation
typedef struct PARAM_T {
    tunable_parameters_t tunable_params;
    int number_fluid_particles_global;
    int number_fluid_particles_local; // Number of non vacant particles not including halo
    int max_fluid_particle_index;     // Max index used in actual particle array
    int number_halo_particles;        // Starting at max_fluid_particle_index
    int max_fluid_particles_local;    // Maximum number of fluid particles
    int number_halo_particles_left;   // Number of halo particles from left neighbor
    int number_halo_particles_right;  // Number of halo particles from right neighbor
    int steps_per_frame;              // Number of simulation steps before updating render node
    float particle_mass; // "mass" of particle so that density is particle count independent
} param_t; // Simulation paramaters

typedef struct render_t {
    float sim_width;
    float sim_height;
    float sim_depth;
    float screen_width;
    float screen_height;
    parameters selected_parameter;
    tunable_parameters_t *node_params; // Holds all nodes paramters including start/end lengths
    tunable_parameters_t *master_params; // Holds parameters shared by all nodes
    int num_compute_procs;
    int num_compute_procs_active; // Number of nodes participating in simulation, user may "remove" nodes at runtime
    bool show_dividers;
    bool pause;
    double last_activity_time; // Used to determine if simulation is being used or not
    bool liquid;
} render_t;

// Struct containing all simulation information
typedef struct FLUID_SIM_T {
    param_t *params;
    AABB_t *water_volume_global;
    AABB_t *boundary_global;
    edge_t *edges;
    oob_t *out_of_bounds;
    neighbor_grid_t *neighbor_grid;      // Neighbor grid setup
    fluid_particles_t *fluid_particles;  // Pointer to fluid_particles SoA
    uint *fluid_particle_indices;        // Index of local fluid particles, used to traverse non vacant particles
    short *fluid_particle_coords;        // (x,y) coordinate array, transfer pixel coords
} fluid_sim_t;

#endif
