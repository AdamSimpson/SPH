#ifndef fluid_geometry_h
#define fluid_geometry_h

typedef struct AABB__ AABB;

#include "fluid.h"

struct AABB__ {
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
}; //Axis aligned bounding box

void constructBoundaryBox(boundary_particle *boundary_particles, AABB* boundary, param *params);

#endif
