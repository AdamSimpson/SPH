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

#ifndef BACKGROUND_GL_H
#define BACKGROUND_GL_H

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

typedef struct background_t {
    GLuint program;

    // Program locations
    GLint position_location;
    GLint tex_coord_location;
    GLint tex_location;

    // Uniforms
    GLuint tex_uniform;

    // Vertex buffer
    GLuint vbo;

    // Element buffer
    GLuint ebo;

    // Pixel dimensions
    int screen_width;
    int screen_height;

    int background_width;
    int background_height;
} background_t;

void init_background(background_t *state, int screen_width, int screen_height);
void create_backround_program(background_t *state);
void create_background_buffers(background_t *state);
void create_background_vertices(background_t *state);
void create_background_texture(background_t *state);
void draw_background(background_t *state);

#endif
