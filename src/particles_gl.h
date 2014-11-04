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

#ifndef PARTICLES_GL_H
#define PARTICLES_GL_H

#include "glfw_utils.h"

typedef struct particles_t {
    // Program handle
    GLuint program;

    // Locations
    GLint position_location;
    GLint color_location;
    GLint diameter_pixels_location;
    GLint radius_world_location;
    GLint view_matrix_location;
    GLint proj_matrix_location;

    // Screen dimensions
    int screen_width;
    int screen_height;

    // buffers
    GLuint vbo;
} particles_t;

// This gets compiled by both c and c++ compilers
#ifdef __cplusplus
extern "C" {
#endif
void init_particles(particles_t *state, int screen_width, int screen_height);
void render_particles(float *points, float diameter_pixels, int num_points, particles_t *state);
void create_particle_shaders(particles_t *state);
void draw_particles(particles_t *state, float diameter_pixels, int num_points);
void create_particle_buffers(particles_t *state);
#ifdef __cplusplus
}
#endif

#endif
