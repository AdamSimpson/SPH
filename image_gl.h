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

#ifndef IMAGE_GL_H
#define IMAGE_GL_H

#include "ogl_utils.h"
#include "stdbool.h"

typedef struct image_t {
    gl_t *gl_state;

    GLuint program;

    // Image file name
    char file_name[50];
    char selected_file_name[50];

    // Program locations
    GLint position_location;
    GLint tex_coord_location;
    GLint tex_location;

    // Uniforms
    GLuint tex_uniform;
    GLuint tex_selected_uniform;

    // Vertex buffer
    GLuint vbo;

    // Element buffer
    GLuint ebo;

    // GL position of lower left corner of image
    float lower_left_x;
    float lower_left_y;

    int image_width;
    int image_height;

    // Determine which image to use
    bool selected;
} image_t;

void init_image(image_t *state, gl_t *gl_state, char *file_name, char *selected_file_name, float lower_left_x, float lower_left_y, int image_width, int image_height);
void create_image_program(image_t *state);
void create_image_buffers(image_t *state);
void create_image_vertices(image_t *state);
void create_image_texture(image_t *state);
void draw_image(image_t *state);
bool inside_image(image_t *state, float gl_x, float gl_y);

#endif
