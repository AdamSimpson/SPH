#ifndef __STRUCTS_H
#define __STRUCTS_H

extern "C"
{
    #include <GL/glew.h>
    #include <GLFW/glfw3.h>
}
#include "stdbool.h"

typedef unsigned int uint;

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

// These parameters are tunable by the render node
struct tunable_parameters_t {
    float rest_density;
    float smoothing_radius;
    float g;
    float k;
    float dq;
    float c;
    float time_step;
    float proc_start;
    float proc_end;
    float mover_center_x;
    float mover_center_y;
    float mover_center_z;
    float mover_radius;
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
