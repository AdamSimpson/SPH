#ifndef fluid_hash_h
#define fluid_hash_h

#include "fluid.h"

//unsigned int hash_val(double x, double y, double z, double h, int hash_size);
uint64_t hash_val(double x, double y, double z, param *params);
void hash_fluid(fluid_particle* fluid_particles, neighbor* neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, param *params);

#endif
