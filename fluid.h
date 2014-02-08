#ifndef fluid_fluid_h
#define fluid_fluid_h

typedef struct FLUID_PARTICLE fluid_particle;
typedef struct NEIGHBOR neighbor;
typedef struct PARAM param;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hash.h"
#include "geometry.h"
#include "communication.h"

// Debug print statement
#define DEBUG 0
#define debug_print(...) \
            do { if (DEBUG) fprintf(stderr, __VA_ARGS__); } while (0)

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct FLUID_PARTICLE {
    float x_prev;
    float y_prev;
    float x;
    float y;
    float v_x;
    float v_y;
    float a_x;
    float a_y;
    float density;
    float density_near;
    float pressure;
    float pressure_near;
    int id; // Id is 'local' index within the fluid particle pointer array
};

struct NEIGHBOR{
    fluid_particle **fluid_neighbors;
    int number_fluid_neighbors;
};

struct PARAM {
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
    int max_bucket_size; // Maximum particles in hash bucket
    int max_neighbors;   // Maximum number of neighbor particles per particle
    int grid_size_x;
    int grid_size_y;
    int number_fluid_particles_global;
    int number_fluid_particles_local; // Number of non vacant particles not including halo
    int max_fluid_particle_index;     // Max index used in actual particle array
    int max_fluid_particles_local;    // Maximum number for max_fluid_particle_index + halo particles
    int number_halo_particles;        // Starting at max_fluid_particle_index
    int length_hash;
}; // Simulation paramaters

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
void collisionImpulse(fluid_particle *p, float norm_x, float norm_y, param *params);
void boundaryConditions(fluid_particle *p, AABB *boundary, param *params);
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                   neighbor *neighbors, n_bucket *hash, AABB* water, int start_x, int number_particles_x, edge *edges, param* params);

void start_simulation();
void calculate_density(fluid_particle *p, fluid_particle *q, float ratio);
void apply_gravity(fluid_particle **fluid_particle_pointers, param *params);
void viscosity_impluses(fluid_particle **fluid_particle_pointers, neighbor* neighbors, param *params);
void predict_positions(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, AABB *boundary_global, param *params);
void double_density_relaxation(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params);
void updateVelocity(fluid_particle *p, param *params);
void updateVelocities(fluid_particle **fluid_particle_pointers, edge *edges, AABB *boundary_global, param *params);
void checkVelocity(float *v_x, float *v_y);


#endif
