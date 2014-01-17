#ifndef MULTI_TEX_H
#define MULTI_TEX_H

#include "GLES2/gl2.h"
#include "egl_utils.h"

#define NUM_TEXTURES 2
#define LEFT 0
#define RIGHT 1

typedef struct
{
    // OpenGL|ES state
    EGL_STATE_T egl_state;

    // Program handle
    GLuint program;

    // Locations
    GLint position_location;
    GLint color_location;

    // buffers
    GLuint vbo;
} STATE_T;

void update_points(float *points, int num_points, STATE_T *state);
void create_shaders(STATE_T *state);
void draw_circles(STATE_T *state, int num_points);

#endif
