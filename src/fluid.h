#ifndef fluid_fluid_h
#define fluid_fluid_h

typedef struct FLUID_PARTICLE fluid_particle_t;
typedef struct NEIGHBOR neighbor_t;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hash.h"
#include "fileio.h"
#include "geometry.h"
#include "communication.h"
#include "sph.h"

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

struct FLUID_PARTICLE {
    double x_star;
    double y_star;
    double z_star;
    double x;
    double y;
    double z;
    double v_x;
    double v_y;
    double v_z;
    double dp_x;
    double dp_y;
    double dp_z;
    double density;
    double lambda;
    int id; // Id is 'local' index within the fluid particle pointer array
};

struct NEIGHBOR {
    unsigned int neighbor_indices[201];
    int number_fluid_neighbors;
};

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
double W(double r, double h);
double del_W(double r, double h);
void XSPH_viscosity(fluid_particle_t *fluid_particles, neighbor_t* neighbors, param_t *params);
void vorticity_confinement(fluid_particle_t *fluid_particles, neighbor_t* neighbors, param_t *params);
void compute_densities(fluid_particle_t *fluid_particles, neighbor_t *neighbors, param_t *params);
void apply_gravity(fluid_particle_t *fluid_particles, param_t *params);
void update_dp_positions(fluid_particle_t *fluid_particles, AABB_t *boundary_global, param_t *params);
void update_positions(fluid_particle_t *fluid_particles, param_t *params);
void calculate_lambda(fluid_particle_t *fluid_particles, neighbor_t *neighbor_grid, param_t *params);
void update_dp(fluid_particle_t *fluid_particles, neighbor_t *neighbor_grid, param_t *params);
void identify_oob_particles(fluid_particle_t *fluid_particles, oob_t *out_of_bounds, param_t *params);
void predict_positions(fluid_particle_t *fluid_particles, AABB_t *boundary_global, param_t *params);
void check_velocity(double *v_x, double *v_y, double *v_z);
void update_velocities(fluid_particle_t *fluid_particles, param_t *params);
void boundary_conditions(fluid_particle_t *fluid_particles, unsigned int i, AABB_t *boudnary);
void initParticles(fluid_particle_t *fluid_particles,
                neighbor_t *neighbors, bucket_t *hash, AABB_t* water, AABB_t* boundary_global,
                int start_x, int number_particles_x, edge_t *edges, param_t* params);


#endif
