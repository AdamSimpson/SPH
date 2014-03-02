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
    char mover_type;

    // Program handle
    GLuint sphere_program;
    GLuint rectangle_program;

    // Locations
    GLint sphere_position_location;
    GLint sphere_center_location;
    GLint sphere_color_location;
    GLint sphere_radius_location;

    // Locations
    GLint rectangle_position_location;
    GLint rectangle_center_location;
    GLint rectangle_color_location;

    // buffers
    GLuint vbo;
} MOVER_T;

void init_mover(MOVER_T *state);
void update_mover(float *center, float *gl_dims, float *color, MOVER_T *state);
void create_sphere_mover_program(MOVER_T *state);
void create_rectangle_mover_program(MOVER_T *state);
void draw_circle_mover(MOVER_T *state, float *center, float radius, float *color);
void draw_rectangle_mover(MOVER_T *state, float *center, float *color);
void create_mover_buffers(MOVER_T *state);

#endif
