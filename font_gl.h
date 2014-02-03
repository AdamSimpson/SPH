#ifndef FONTS_GL_H
#define FONTS_GL_H

#include "renderer.h"

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

typedef struct
{
    struct FONScontext* fs;
    int font_normal;

    int screen_width;
    int screen_height;
} FONT_T;

void init_font(FONT_T *font_state, int screen_width, int screen_height);
void render_parameters(FONT_T *font_state, parameters selected_param, double gravity, double viscosity, double density, double pressure, double elasticity);
void render_fps(FONT_T *font_state, double fps);
void remove_font(FONT_T *font_state);


#endif
