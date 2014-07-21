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

#ifndef FONTS_GL_H
#define FONTS_GL_H

#include "renderer.h"

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

#include <ft2build.h>
#include FT_FREETYPE_H

typedef struct text_vert_t {
    GLfloat x;
    GLfloat y;
    GLfloat s;
    GLfloat t;
    GLfloat r;
    GLfloat g;
    GLfloat b;      
} text_vert_t;

// Structure to hold cache glyph information
typedef struct char_info_t {
    float ax; // advance.x
    float ay; // advance.y

    float bw; // bitmap.width
    float bh; // bitmap.height

    float bl; // bitmap left
    float bt; // bitmap top

    float tx; // x offset of glyph in texture coordinates
    float ty; // y offset of glyph in texture coordinates
} char_info_t;

typedef struct font_t
{
    FT_Library ft;
    FT_Face face;

    // Font shader program
    GLuint program;

    // Program locations
    GLint coord_location;
    GLint tex_location;
    GLint color_location;
    
    // Uniforms
    GLuint tex_uniform;

    // VBO
    GLuint vbo;

    int screen_width;
    int screen_height;
    
    // Font atlas
    char_info_t char_info[128];
    int atlas_width;
    int atlas_height;
} font_t;

void create_font_program(font_t *state);
void create_font_buffers(font_t *state);
void create_font_atlas(font_t *state);
void init_font(font_t *state, int screen_width, int screen_height);
void render_fps(font_t *state, float fps);
void render_parameters(font_t *state, render_t *render_state);
int add_text_coords(font_t *state, char *text, text_vert_t* verts, float *color, float x, float y, float sx, float sy);
void render_all_text(font_t *state, render_t *render_state, double fps);

#endif
