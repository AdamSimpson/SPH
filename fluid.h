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
#include "fileio.h"
#include "geometry.h"
#include "communication.h"

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct FLUID_PARTICLE {
    double x_prev;
    double y_prev;
    double x;
    double y;
    double v_x;
    double v_y;
    double a_x;
    double a_y;
    double density;
    double density_near;
    double pressure;
    double pressure_near;
    int id; // Id is 'local' index within the fluid particle pointer array
};

struct NEIGHBOR{
    fluid_particle* fluid_neighbors[300];
    int number_fluid_neighbors;
};

struct PARAM {
    double rest_density;
    double mass_particle;
    double spacing_particle;
    double smoothing_radius;
    double g;
    double time_step;
    double alpha;
    double surface_tension;
    double speed_sound;
    double node_start_x; // left x position of node partition
    double node_end_x;   // right x position of node partition
    int grid_size_x;
    int grid_size_y;
    int number_fluid_particles_global;
    int number_fluid_particles_local; // Number of non vacant particles
    int max_fluid_particle_index;     // Max index used in actual particle array
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
double lap_W_visc(const double r, const double h);
double W_dens(const double r, const double h);
double del_W_pressure(const double r, const double h);
void collisionImpulse(fluid_particle *p, int norm_x, int norm_y);
void boundaryConditions(fluid_particle *p, AABB *boundary, param *params);
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                   neighbor *neighbors, n_bucket *hash, AABB* water, int start_x, int number_particles_x, edge *edges, param* params);

void calculate_density(fluid_particle *p, fluid_particle *q, param *params);
void apply_gravity(fluid_particle **fluid_particle_pointers, param *params);
void viscosity_impluses(fluid_particle **fluid_particle_pointers, neighbor* neighbors, param *params);
void predict_positions(fluid_particle **fluid_particle_pointers, AABB *boundary_global, param *params);
void double_density_relaxation(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params);
void updateVelocity(fluid_particle *p, param *params);
void updateVelocities(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, edge *edges, AABB *boundary_global, param *params);

#endif
