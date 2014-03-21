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
