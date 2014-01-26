#ifndef fluid_geometry_h
#define fluid_geometry_h

typedef struct AABB__ AABB;

#include "fluid.h"
#include "communication.h"

struct AABB__ {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
}; //Axis aligned bounding box

double min(double a, double b);
double max(double a, double b);
int sgn(double x);
void constructFluidVolume(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, AABB* fluid, int start_x, int number_particles_x, edge *edges, param *params);
void setParticleNumbers(AABB *boundary_global, AABB *fluid_global, edge *edges, oob *out_of_bounds, int number_particles_x, param *params);
void partitionProblem(AABB *boundary_global, AABB *fluid_global, int *x_start, int *number_particles_x, param *params);
void checkPartition(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, double partition_time, param *params);
#endif
