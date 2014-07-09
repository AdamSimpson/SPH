#ifndef FLUID_CUDA_H
#define FLUID_CUDA_H

#include "fluid.h"

void double_density_relaxation_gpu(fluid_particle **fluid_particle_pointers, uint *particle_ids, uint *start_indexes, uint *end_indexes, param *params);
void calculate_pressures_gpu(fluid_particle **fluid_particle_pointers, param *params);
void updateVelocities_gpu(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params);
void predict_positions_gpu(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params);
void hash_particles_gpu(fluid_particle **fluid_particle_pointers, uint *hash_values, uint *particle_ids, uint *starts, uint *ends, param *params);
void calculate_density_gpu(fluid_particle **fluid_particle_pointers, uint *start_indexes, uint *end_indexes, uint *particle_ids, param *params);
void viscosity_impluses_gpu(fluid_particle **fluid_particle_pointers, uint *particle_ids, uint *start_indexes, uint *end_indexes, param *params);
void sort_hash_gpu(uint *d_particle_ids, uint *d_hash_values, param *params);
void apply_gravity_gpu(fluid_particle **fluid_particle_pointers, param *params);

#endif
