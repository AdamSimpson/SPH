#ifndef fluid_geometry_h
#define fluid_geometry_h

typedef struct AABB AABB_t;

#include "fluid.h"
#include "communication.h"
#include "sph.h"

struct AABB {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
};//Axis aligned bounding box

double min(double a, double b);
double max(double a, double b);
int sgn(double x);
void constructFluidVolume(fluid_particle_t *fluid_particles, AABB_t* fluid, int start_x, int number_particles_x, edge_t *edges, param_t *params);
void setParticleNumbers(AABB_t *boundary_global, AABB_t *fluid_global, edge_t *edges, oob_t *out_of_bounds, int number_particles_x, param_t *params);
void partitionProblem(AABB_t *boundary_global, AABB_t *fluid_global, int *x_start, int *number_particles_x, param_t *params);
void checkPartition(fluid_particle_t *fluid_particles, oob_t *out_of_bounds, param_t *params);
#endif
