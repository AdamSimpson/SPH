#ifndef SPH_H
#define SPH_H

typedef struct PARAM param_t;

#include "fluid.h"
#include "hash.h"
#include "geometry.h"
#include "fileio.h"
#include "communication.h"

struct PARAM {
    double rest_density;
    double smoothing_radius;
    double g;
    double time_step;
    double k;
    double dq;
    double c;
    double node_start_x; // left x position of node partition
    double node_end_x;   // right x position of node partition
    int grid_size_x;
    int grid_size_y;
    int grid_size_z;
    int number_fluid_particles_global;
    int number_fluid_particles_local;  // Number of non vacant particles
    int max_fluid_particles_local;     // Maximum number for max_fluid_particle_index + halo particles
    int number_halo_particles_left;    // Starting at max_fluid_particle_index
    int number_halo_particles_right;
    int number_steps;
    int length_hash;
    int rank;
    int nprocs;
}; // Simulation paramaters

#endif
