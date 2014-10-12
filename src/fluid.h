/*
The MIT License (MIT)

Copyright (c) 2014 Adam Simpson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#ifndef fluid_fluid_h
#define fluid_fluid_h

typedef struct FLUID_PARTICLE_T fluid_particle_t;
typedef struct PARAM_T param_t;
typedef struct TUNABLE_PARAMETERS_T tunable_parameters_t;
typedef struct FLUID_SIM_T fluid_sim_t;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hash.h"
#include "setup.h"
#include "communication.h"

// Debug print statement
#define DEBUG 0
#define debug_print(...) \
            do { if (DEBUG) fprintf(stderr, __VA_ARGS__); } while (0)


////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

// Standard fluid particle paramaters
struct FLUID_PARTICLE_T {
    float x_star;
    float y_star;
    float x;
    float y;
    float v_x;
    float v_y;
    float density;
    float dp_x;
    float dp_y;
    float lambda;
    int id; // Id is 'local' index within the fluid particle pointer array
};

// These parameters are tunable by the render node
struct TUNABLE_PARAMETERS_T {
    float rest_density;
    float smoothing_radius;
    float g;
    float k;
    float k_near;
    float k_spring;
    float sigma;
    float beta;
    float time_step;
    float node_start_x;
    float node_end_x;
    float mover_center_x;
    float mover_center_y;
    float mover_width;
    float mover_height;
    char mover_type;
    char kill_sim;
    char active;
};

// Full parameters struct for simulation
struct PARAM_T {
    tunable_parameters_t tunable_params;
    int number_fluid_particles_global;
    int number_fluid_particles_local; // Number of non vacant particles not including halo
    int max_fluid_particle_index;     // Max index used in actual particle array
    int number_halo_particles;        // Starting at max_fluid_particle_index
    int max_fluid_particles_local;    // Maximum number of fluid particles
    int number_halo_particles_left;   // Number of halo particles from left neighbor
    int number_halo_particles_right;  // Number of halo particles from right neighbor
    int steps_per_frame;              // Number of simulation steps before updating render node
    float particle_mass; // "mass" of particle so that density is particle count independent
}; // Simulation paramaters

// Struct containing all simulation information
struct FLUID_SIM_T {
    param_t *params;
    AABB_t *water_volume_global;
    AABB_t *boundary_global;
    edge_t *edges;
    oob_t *out_of_bounds;
    neighbor_grid_t *neighbor_grid;  // Neighbor grid setup
    fluid_particle_t *fluid_particles; // Storage for fluid particles
    short *fluid_particle_coords;    // (x,y) coordinate array, transfer pixel coords
    fluid_particle_t **fluid_particle_pointers;  //pointer array used to traverse non vacant particles
};

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
//void collisionImpulse(fluid_particle *p, float norm_x, float norm_y, param *params);
void boundaryConditions(fluid_particle_t *p, fluid_sim_t *fluid_sim);
void init_sim_particles(fluid_sim_t *fluid_sim, float start_x, int number_particles_x);
void start_simulation();
void calculate_density(fluid_particle_t *p, fluid_particle_t *q, float h, float mass);
void apply_gravity(fluid_sim_t *fluid_sim);
void updateVelocity(fluid_particle_t *p, param_t *params);
void updateVelocities(fluid_sim_t *fluid_sim);
void checkVelocity(float *v_x, float *v_y);
void identify_oob_particles(fluid_sim_t *fluid_sim);
float del_W(float r, float h);
float W(float r, float h);
void predict_positions(fluid_sim_t *fluid_sim);
void update_dp_positions(fluid_sim_t *fluid_sim);
void update_positions(fluid_sim_t *fluid_sim);
void calculate_lambda(fluid_sim_t *fluid_sim);
void update_dp(fluid_sim_t *fluid_sim);
void compute_densities(fluid_sim_t *fluid_sim);
void XSPH_viscosity(fluid_sim_t *fluid_sim);
void vorticity_confinement(fluid_sim_t *fluid_sim);

#endif
