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

typedef struct FLUID_PARTICLE fluid_particle;
typedef struct NEIGHBOR neighbor;
typedef struct PARAM param;
typedef struct TUNABLE_PARAMETERS tunable_parameters;

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "geometry.h"
#include "communication.h"

// Debug print statement
#define DEBUG 0
#define debug_print(...) \
            do { if (DEBUG) fprintf(stderr, __VA_ARGS__); } while (0)

// MPI doesn't have a C enum type
// Defines will be ok for our use
#define SPHERE_MOVER 0
#define RECTANGLE_MOVER 1

////////////////////////////////////////////////
// Structures
////////////////////////////////////////////////

// Standard fluid particle paramaters
struct FLUID_PARTICLE {
    float x_prev;
    float y_prev;
    float x;
    float y;
    float v_x;
    float v_y;
    float a_x;
    float a_y;
    float density;
    float density_near;
    float pressure;
    float pressure_near;
    int id; // Id is 'local' index within the fluid particle pointer array
};

struct NEIGHBOR{
    fluid_particle **fluid_neighbors;
    int number_fluid_neighbors;
};

// These parameters are tunable by the render node
struct TUNABLE_PARAMETERS {
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
struct PARAM {
    tunable_parameters tunable_params;
    int number_fluid_particles_global;
    int number_fluid_particles_local; // Number of non vacant particles not including halo
    int max_fluid_particle_index;     // Max index used in actual particle array
    int number_halo_particles;        // Starting at max_fluid_particle_index
    float grid_spacing;
    float grid_size_x;
    float grid_size_y;
}; // Simulation paramaters

////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
//void collisionImpulse(fluid_particle *p, float norm_x, float norm_y, param *params);
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                   AABB_t *water, int start_x, int number_particles_x, 
		   edge_t *edges, int max_fluid_particles_local, float spacing, param* params);
void start_simulation();
void identify_oob_particles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, oob_t *out_of_bounds, AABB_t *boundary_global, param *params);

#endif
