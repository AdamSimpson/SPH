#ifndef FONTS_GL_H
#define FONTS_GL_H

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

typedef struct
{
    struct FONScontext* fs;
    int font_normal;

} FONT_T;

void init_font(FONT_T *font_state);
void render_font(FONT_T *font_state);
void remove_font(FONT_T *font_state);


#endif
