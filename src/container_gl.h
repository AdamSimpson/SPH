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

#ifndef CONTAINER_GL_H
#define CONTAINER_GL_H

#include "glfw_utils.h"

typedef struct container_t {
    GLuint program;

    // Program locations
    GLint position_location;
    GLint normal_location;
    GLint tex_coord_location;

    // Uniforms
    GLuint color_location;

    GLint global_matrix_index;

    // Vertex objects
    GLuint vbo;
    GLuint vao;

    // Pixel dimensions
    int screen_width;
    int screen_height;

    float width;
    float height;
    float depth;

    // Uniform objects
    GLuint ubo;
} container_t;

// This gets compiled by both c and c++ compilers
#ifdef __cplusplus
extern "C" {
#endif
void init_container(container_t *state, int screen_width, int screen_height);
void create_container_program(container_t *state);
void create_container_buffers(container_t *state);
void create_container_vertices(container_t *state);
void render_container(container_t *state);
#ifdef __cplusplus
}
#endif

#endif
