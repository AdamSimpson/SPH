#ifndef fluid_hash_h
#define fluid_hash_h

#include <stdbool.h>

typedef struct N_BUCKET n_bucket;

#include "fluid.h"

struct N_BUCKET {
    fluid_particle *fluid_particles[201];
    unsigned int number_fluid;
    bool hashed;
}; // neighbor 'bucket' for hash value

unsigned int hash_val(double x, double y, double z, double h, int hash_size);
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket * hash, param *params);
void hash_halo(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket *hash, param *params);

#endif

