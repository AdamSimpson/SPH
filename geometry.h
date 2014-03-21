#ifndef fluid_geometry_h
#define fluid_geometry_h

typedef struct AABB_T AABB_t;

#include "fluid.h"
#include "communication.h"

struct AABB_T {
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
void partitionProblem(AABB_t *boundary_global, AABB_t *fluid_global, int *x_start, int *length_x, float spacing, param *params);
void setParticleNumbers(AABB_t *boundary_global, AABB_t *fluid_global, edge_t *edges, oob_t *out_of_bounds, int number_particles_x, float spacing, param *params);

void constructFluidVolume(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, AABB_t* fluid, int start_x, 
                          int number_particles_x, edge_t *edges, float spacing, param *params);

#endif
