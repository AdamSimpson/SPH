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

#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "hash_sort.h"
#include "setup.h"
#include "communication.h"

// Debug print statement
#define DEBUG 0
#define debug_print(...) \
            do { if (DEBUG) fprintf(stderr, __VA_ARGS__); } while (0)


////////////////////////////////////////////////
// Function prototypes
////////////////////////////////////////////////
//void collisionImpulse(fluid_particle *p, float norm_x, float norm_y, param *params);
void boundary_conditions(uint p_index, fluid_sim_t *fluid_sim);
void start_simulation();
void calculate_density(uint p_index, uint q_index, float h, float mass);
void apply_gravity(fluid_sim_t *fluid_sim);
void update_velocities(fluid_sim_t *fluid_sim);
void check_velocity(float *v_x, float *v_y, float *v_z);
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
