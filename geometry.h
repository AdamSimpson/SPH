#ifndef fluid_geometry_h
#define fluid_geometry_h

typedef struct AABB__ AABB;

#include "fluid.h"
#include "communication.h"

struct AABB__ {
    float min_x;
    float max_x;
    float min_y;
    float max_y;
    float min_z;
    float max_z;
}; //Axis aligned bounding box

float min(float a, float b);
float max(float a, float b);
int sgn(float x);
void constructFluidVolume(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, AABB* fluid, int start_x, int number_particles_x, edge *edges, param *params);
void setParticleNumbers(AABB *boundary_global, AABB *fluid_global, edge *edges, oob *out_of_bounds, int number_particles_x, param *params);
void partitionProblem(AABB *boundary_global, AABB *fluid_global, int *x_start, int *number_particles_x, param *params);
void checkPartition(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, double partition_time, param *params);
#endif
