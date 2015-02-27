#ifndef fluid_hash_h
#define fluid_hash_h

typedef struct BUCKET bucket_t;

#include <stdbool.h>
#include "fluid.h"
#include "geometry.h"
#include "sph.h"

struct BUCKET {
    fluid_particle_t *fluid_particles[201];
    unsigned int number_fluid;
    bool hashed;
}; // neighbor 'bucket' for hash value

unsigned int hash_val(double x, double y, double z, param_t *params);
void hash_fluid(fluid_particle_t *fluid_particles, neighbor_t *neighbors, bucket_t *hash, AABB_t *boundary, param_t *params);
void hash_halo(fluid_particle_t *fluid_particles, neighbor_t *neighbors, bucket_t *hash, AABB_t *boundary, param_t *params);

#endif
