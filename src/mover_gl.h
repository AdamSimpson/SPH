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

#ifndef MOVER_GL_H
#define MOVER_GL_H

#ifdef GLFW
    #include "glfw_utils.h"
#else
    #include "GLES2/gl2.h"
    #include "egl_utils.h"
#endif

typedef struct mover_t
{
    // Program handle
    GLuint sphere_program;
    GLuint rectangle_program;

    // Locations
    GLint sphere_position_location;
    GLint sphere_center_location;
    GLint sphere_tex_coord_location;
    GLint sphere_color_location;
    GLint sphere_radius_location;

    // buffers
    GLuint vbo;
} mover_t;

void init_mover(mover_t *state);
void render_mover(float *center, float *gl_dims, float *color, mover_t *state);
void create_sphere_mover_program(mover_t *state);
void draw_circle_mover(mover_t *state, float *center, float radius, float *color);
void create_mover_buffers(mover_t *state);

#endif
