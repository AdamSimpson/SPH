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
    FT_GlyphSlot g;

    // Font shader program
    GLuint program;

    // Coordinate location
    GLint coord_location;

    // texture uniform
    GLuint tex_uniform;

    // VBO
    GLuint vbo;

    int screen_width;
    int screen_height;
} FONT_T;

void create_shaders(FONT_T *state);
void create_buffers(FONT_T *state);


#endif
