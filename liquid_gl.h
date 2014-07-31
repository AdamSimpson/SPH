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

#ifndef LIQUID_GL_H
#define LIQUID_GL_H

#ifdef GLFW
    #include "glfw_utils.h"
#else
    #include "GLES2/gl2.h"
    #include "egl_utils.h"
#endif

typedef struct liquid_t {
    // Program handle
    GLuint program;

    // Program to render texture to screen
    GLuint tex_program;

    // Program for vert gaussian blur
    GLuint vert_blur_program;

    // Program for horizontal Gaussian blur
    GLuint horz_blur_program;

    // Render to low rez tex Locations
    GLint position_location;
    GLint diameter_pixels_location;

    // Gaussian locations
    GLint vert_blur_position_location;
    GLint vert_blur_tex_coord_location;
    GLint vert_blur_tex_location;
    GLint horz_blur_position_location;
    GLint horz_blur_tex_coord_location;
    GLint horz_blur_tex_location;

    // Render tex to quad Locations
    GLint tex_position_location;
    GLint tex_location;
    GLint tex_coord_location;

    // Uniforms
    GLuint tex_uniform;
    GLuint blur_horz_tex_uniform;

    // Screen dimensions
    int screen_width;
    int screen_height;

    // buffers
    GLuint vbo;
    GLuint tex_vbo;
    GLuint tex_ebo;

    // buffers for meta ball render to texture
    GLuint frame_buffer;
    GLuint tex_color_buffer;
    GLuint blur_horz_color_buffer;

    GLuint reduction;
} liquid_t;

void init_liquid(liquid_t *state, int screen_width, int screen_height);
void render_liquid(float *points, float diameter_pixels, int num_points, liquid_t *state);
void create_liquid_shaders(liquid_t *state);
void draw_liquid(liquid_t *state, float diameter_pixels, int num_points);
void create_liquid_buffers(liquid_t *state);
void create_texture_verticies(liquid_t *state);

#endif
