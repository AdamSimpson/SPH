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

#ifndef CURSOR_GL_H
#define CURSOR_GL_H

#include "ogl_utils.h"
#include "stdbool.h"

typedef struct cursor_t {
    gl_t *gl_state;

    GLuint program;

    char file_name[40];

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

    // Position of cursor
    float center_x;
    float center_y;

    // Pixel dimensions of cursor image
    int cursor_width;
    int cursor_height;
} cursor_t;

void init_cursor(cursor_t *state, gl_t *gl_state, char *file_name, int cursor_width, int cursor_height);
void create_cursor_program(cursor_t *state);
void create_cursor_buffers(cursor_t *state);
void set_cursor_position(cursor_t *state, float gl_x, float gl_y);
void create_cursor_texture(cursor_t *state);
void draw_cursor(cursor_t *state);

#endif
