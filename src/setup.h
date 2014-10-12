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

#ifndef fluid_setup_h
#define fluid_setup_h

typedef struct AABB_T AABB_t;

#include "fluid.h"
#include "communication.h"

struct AABB_T {
    float min_x;
    float max_x;
    float min_y;
    float max_y;
    float min_z;
    float max_z;
}; //Axis aligned bounding box

float min(float a, float b);
float max(float a, float b);
int sgn(float x);
void set_particle_numbers(fluid_sim_t *fluid_sim, int number_particles_x, float spacing);
void partition_geometry(fluid_sim_t *fluid_sim, float *x_start, int *number_particles_x, float spacing);
void partition_simulation(fluid_sim_t *fluid_sim, float *x_start, int *number_particles_x);
void construct_fluid_volume(fluid_sim_t *fluid_sim, float start_x, int number_particles_x);
void alloc_sim_structs(fluid_sim_t *fluid_sim);
void free_sim_structs(fluid_sim_t *fluid_sim);
void init_params(fluid_sim_t *fluid_sim);
void sync_initial_params(fluid_sim_t *fluid_sim);
void alloc_sim(fluid_sim_t *fluid_sim);
void free_sim_memory(fluid_sim_t *fluid_sim);

#endif
