#ifndef fluid_hash_h
#define fluid_hash_h

#include "fluid.h"

//unsigned int hash_val(double x, double y, double z, double h, int hash_size);
uint64_t hash_val(double x, double y, double z, param *params);
void hash_boundary(boundary_particle *boundary_particles, uint2 *hash, uint2 *boundary_hash_positions,  param *params);
void hash_fluid(fluid_particle* fluid_particles, boundary_particle *boundary_particles, neighbor* neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, uint2 *boundary_hash, uint2 *boundary_hash_positions, param *params);

#endif
