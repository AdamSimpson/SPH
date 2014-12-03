#ifndef __STRUCTS_H
#define __STRUCTS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "stdbool.h"

typedef unsigned int uint;

// Poor choices have been made
#ifndef __cplusplus
typedef struct gl_t gl_t;
typedef struct world_t world_t;
typedef struct fluid_particles_t fluid_particles_t;
typedef struct edge_t edge_t; 
typedef struct oob_t oob_t;
typedef struct AABB_t AABB_t;
typedef struct neighbor_t neighbor_t;
typedef struct neighbor_grid_t neighbor_grid_t;
typedef struct tunable_parameters_t tunable_parameters_t;
typedef struct param_t param_t;
typedef struct render_t render_t;
typedef struct fluid_sim_t fluid_sim_t;
#endif

// gl_t must be lowercase!!! Something else is using GL_T (?)
struct gl_t {
    int screen_width;
    int screen_height;

    GLFWwindow* window;
};

struct world_t {
    // Screen dimensions
    int screen_width;
    int screen_height;

    float eye_position[3];
    float look_at[3];
};

// Standard fluid particle paramaters
struct fluid_particles_t {
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
};

// Particles that are within 2*h distance of node edge
struct edge_t {
    int max_edge_particles;
    uint *edge_indices_left;
    uint *edge_indices_right;
    int number_edge_particles_left;
    int number_edge_particles_right;
};

// Particles that have left the node
struct oob_t {
    int max_oob_particles;
    uint *oob_index_indices_left; // Indicies in particle index array for particles traveling left
    uint *oob_index_indices_right;
    int number_oob_particles_left; // Number of OOB sending left
    int number_oob_particles_right;
    uint *vacant_indices; // Indicies in global particle array that are vacant
    int number_vacancies;
};

struct AABB_t {
    float min_x;
    float max_x;
    float min_y;
    float max_y;
    float min_z;
    float max_z;
}; //Axis aligned bounding box

// Bucket to hold each particles nearest neighbors
struct neighbor_t{
    uint *fluid_neighbors; // Index in global particle array of neighbor
    int number_fluid_neighbors;
};

struct neighbor_grid_t {
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
};

// enum of displayed parameter values
typedef enum {
    MIN = 0,
    GRAVITY = MIN,
    SMOOTH,
    DENSITY,
    K,
    DQ,
    VISCOSITY,
    MAX = VISCOSITY
} parameters;

// These parameters are tunable by the render node
struct tunable_parameters_t {
    float rest_density;
    float smoothing_radius;
    float g;
    float k;
    float dq;
    float c;
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
};

// Full parameters struct for simulation
struct param_t {
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
}; // Simulation paramaters

struct render_t {
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
    bool view_controls; // When shift is held mouse controls view
    bool pause;
    double last_activity_time; // Used to determine if simulation is being used or not
    world_t *world;
};

// Struct containing all simulation information
struct fluid_sim_t {
    param_t *params;
    AABB_t *water_volume_global;
    AABB_t *boundary_global;
    edge_t *edges;
    oob_t *out_of_bounds;
    neighbor_grid_t *neighbor_grid;      // Neighbor grid setup
    fluid_particles_t *fluid_particles;  // Pointer to fluid_particles SoA
    uint *fluid_particle_indices;        // Index of local fluid particles, used to traverse non vacant particles
    short *fluid_particle_coords;        // (x,y) coordinate array, transfer pixel coords
};

#endif
