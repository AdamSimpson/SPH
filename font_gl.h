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

typedef struct
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
    GLuint color_uniform;

    // VBO
    GLuint vbo;

    int screen_width;
    int screen_height;
} FONT_T;

void create_font_program(FONT_T *state);
void create_font_buffers(FONT_T *state);
void create_font_textures(FONT_T *state);
void init_font(FONT_T *state, int screen_width, int screen_height);
void render_fps(FONT_T *state, double fps);

#endif
