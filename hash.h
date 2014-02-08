#ifndef fluid_hash_h
#define fluid_hash_h

#include <stdbool.h>

typedef struct BUCKET_T bucket_t;
typedef struct NEIGHBOR_GRID_T neighbor_grid_t;

#include "fluid.h"

struct BUCKET_T {
    fluid_particle **fluid_particles;
    unsigned int number_fluid;
}; // neighbor 'bucket' for hash value

struct NEIGHBOR_GRID_T {
    float spacing;  // Spacing between buckets
    unsigned int size_x; // Number of buckets in x
    unsigned int size_y; // Number of buckets in y
    neighbor *neighbors; // Particle neighbor buckets
    bucket_t *grid_buckets; // Grid to place hashed particles into
    unsigned int max_neighbors; // Maximum neighbors allowed for each particle
    unsigned int max_bucket_size; // Maximum particles in hash bucket
};

unsigned int hash_val(float x, float y, neighbor_grid_t *grid, param *params);
void hash_fluid(fluid_particle **fluid_particle_pointers, neighbor_grid_t *grid, param *params, bool compute_density);
void hash_halo(fluid_particle **fluid_particle_pointers,  neighbor_grid_t *grid, param *params, bool compute_density);

#endif

