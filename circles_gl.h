#ifndef CIRCLES_GL_H
#define CIRCLES_GL_H

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
    GLint color_location;
    GLint radius_location;

    // buffers
    GLuint vbo;
} CIRCLE_T;

void init_circles(CIRCLE_T *state);
void update_points(float *points, int num_points, CIRCLE_T *state);
void update_mover_point(float *point, float radius, CIRCLE_T *state);
void create_shaders(CIRCLE_T *state);
void draw_circles(CIRCLE_T *state, int num_points);
void draw_circle_mover(CIRCLE_T *state, float radius);
void create_buffers(CIRCLE_T *state);

#endif
