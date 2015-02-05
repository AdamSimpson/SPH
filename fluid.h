#ifndef fluid_fluid_h
#define fluid_fluid_h

typedef struct FLUID_PARTICLE fluid_particle_t;
typedef struct NEIGHBOR neighbor_t;
typedef struct PARAM param_t;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hash.h"
#include "fileio.h"
#include "geometry.h"
#include "communication.h"

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct FLUID_PARTICLE {
    float x_star;
    float y_star;
    float z_star;
    float x;
    float y;
    float z;
    float v_x;
    float v_y;
    float v_z;
    float dp_x;
    float dp_y;
    float dp_z;
    float density;
    float lambda;
    int id; // Id is 'local' index within the fluid particle pointer array
};

struct NEIGHBOR {
    fluid_particle_t* fluid_neighbors[201];
    int number_fluid_neighbors;
};

struct PARAM {
    double rest_density;
    double spacing_particle;
    double smoothing_radius;
    double g;
    double time_step;
    double k;
    double dq;
    double c;
    double node_start_x; // left x position of node partition
    double node_end_x;   // right x position of node partition
    int grid_size_x;
    int grid_size_y;
    int grid_size_z;
    int number_fluid_particles_global;
    int number_fluid_particles_local; // Number of non vacant particles
    int max_fluid_particles_local;    // Maximum number for max_fluid_particle_index + halo particles
    int number_halo_particles;        // Starting at max_fluid_particle_index
    int max_node_difference;
    int number_steps;
    int length_hash;
    int rank;
    int nprocs;
}; // Simulation paramaters

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
float W(float r, float h);
float del_W(float r, float h);
void XSPH_viscosity(fluid_particle_t *fluid_particles, neighbor_t* neighbors, param_t *params);
void compute_densities(fluid_particle_t *fluid_particles, param_t *params);
void apply_gravity(fluid_particle_t *fluid_particles, param_t *params);
void update_dp_positions(fluid_particle_t *fluid_particles, param_t *params);
void update_positions(fluid_particle_t *fluid_particles, param_t *params);
void calculate_lambda(fluid_particle_t *fluid_particles, neighbor_t *neighbor_grid, param_t *params);
void update_dp(fluid_particle_t *fluid_particles, neighbor_t *neighbor_grid, param_t *params);
void identify_oob_particles(fluid_particle_t *fluid_particles, oob_t *out_of_bounds, param_t *params);
void predict_positions(fluid_particle_t *fluid_particles, param_t *params);
void check_velocity(float *v_x, float *v_y, float *v_z);
void update_velocities(fluid_particle_t *fluid_particles, param_t *params);
void boundary_conditions(fluid_particle_t *fluid_particles, unsigned int i, AABB_t *boudnary);
void initParticles(fluid_particle_t *fluid_particles,
                neighbor_t *neighbors, bucket_t *hash, AABB_t* water,
                int start_x, int number_particles_x, edge_t *edges, param_t* params);


#endif
