#ifndef fluid_fluid_h
#define fluid_fluid_h

typedef struct BOUNDARY_PARTICLE boundary_particle;
typedef struct FLUID_PARTICLE fluid_particle;
typedef struct NEIGHBOR neighbor;
typedef struct PARAM param;
typedef struct UINT2 uint2;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include "hash.h"
#include "fileio.h"
#include "geometry.h"

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct UINT2 {
    uint64_t x;
    uint64_t y;
};

struct BOUNDARY_PARTICLE {
    double x;
    double y;
    double z;
    double n_x;
    double n_y;
    double n_z;
};

struct FLUID_PARTICLE {
    double density;
    double pressure;
    double x;
    double y;
    double z;
    double v_x;
    double v_y;
    double v_z;
    double v_half_x;
    double v_half_y;
    double v_half_z;
    double a_x;
    double a_y;
    double a_z;
};

struct NEIGHBOR{
    fluid_particle* fluid_neighbors[65];
    boundary_particle* boundary_neighbors[35];
    int number_fluid_neighbors;
    int number_boundary_neighbors;
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
    int number_particles;
    int number_fluid_particles;
    int number_boundary_particles;
    int number_steps;
    int length_fluid_hash;
    int length_boundary_hash;
    int grid_size_x;
    int grid_size_y;
    int grid_size_z;
}; // Simulation paramaters

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////

double W(fluid_particle *p, fluid_particle *q, double h);
double del_W(fluid_particle *p, fluid_particle *q, double h);
double computeDensity(fluid_particle *p, fluid_particle *q, param *params);
double computePressure(fluid_particle *p, param *params);
void updatePressures(fluid_particle *fluid_particles, neighbor *neighbors, param *params);
void computeAcceleration(fluid_particle *p, fluid_particle *q, param *params);
void computeBoundaryAcceleration(fluid_particle *p, boundary_particle *q, param *params);
void updateAccelerations(fluid_particle *fluid_particles, neighbor *neighbors, param *params);
void updateParticle(fluid_particle *p, param *params);
void updatePositions(fluid_particle *fluid_particles, param *params);
void eulerStart(fluid_particle* fluid_particles, boundary_particle *boundary_particles, neighbor *neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, uint2 *boundary_hash, uint2 *boundary_hash_positions, param *params);
void initParticles(fluid_particle* fluid_particles, boundary_particle* boundary_particles, neighbor *neighbors, uint2 *hash, uint2 *fluid_hash_positions, uint2 *boundary_hash, uint2 *boundary_hash_positions, AABB* water, AABB* boundary, param* params);
#endif
