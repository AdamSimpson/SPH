#ifndef fluid_hash_h
#define fluid_hash_h

#include <stdbool.h>

typedef struct N_BUCKET n_bucket;

#include "fluid.h"

struct N_BUCKET {
    fluid_particle **fluid_particles;
    unsigned int number_fluid;
}; // neighbor 'bucket' for hash value

unsigned int hash_val(double x, double y, param *params);
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket * hash, param *params, bool compute_density);
void hash_halo(fluid_particle **fluid_particle_pointers, neighbor *neighbors, n_bucket *hash, param *params, bool compute_density);

#endif

