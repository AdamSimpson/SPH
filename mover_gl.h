#ifndef MOVER_GL_H
#define MOVER_GL_H

#ifdef GLFW
    #include "glfw_utils.h"
#else
    #include "GLES2/gl2.h"
    #include "egl_utils.h"
#endif

typedef struct
{
    // Program handle
    GLuint program;

    // Locations
    GLint position_location;

    GLint center_location;
    GLint color_location;
    GLint radius_location;

    // buffers
    GLuint vbo;
} MOVER_T;

void init_mover(MOVER_T *state);
void update_mover(float *center, float radius, float *gl_dims, float *color, MOVER_T *state);
void create_mover_shaders(MOVER_T *state);
void draw_circle_mover(MOVER_T *state, float *center, float radius, float *color);
void create_mover_buffers(MOVER_T *state);

#endif
