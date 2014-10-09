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

#ifndef RECTANGLE_GL_H
#define RECTANGLE_GL_H

#include "glfw_utils.h"

typedef struct rectangle_t
{
    // Program handle
    GLuint rectangle_program;

    // Locations
    GLint color_location;
    GLint position_location;

    // buffers
    GLuint vbo;
} rectangle_t;

void init_rectangle(rectangle_t *state);
void render_rectangle(rectangle_t *state, float *center, float *gl_dims, float *color);
void create_rectangle_program(rectangle_t *state);
void draw_rectangle(rectangle_t *state, float *center, float *color);

void create_rectangle_buffers(rectangle_t *state);

#endif
