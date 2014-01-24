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
    // OpenGL state
    GL_STATE_T gl_state;

    // Program handle
    GLuint program;

    // Locations
    GLint position_location;
    GLint color_location;

    // buffers
    GLuint vbo;
} STATE_T;

inline void check();
void showlog(GLint shader);
void update_points(float *points, int num_points, STATE_T *state);
void create_shaders(STATE_T *state);
void draw_circles(STATE_T *state, int num_points);
void create_buffers(STATE_T *state);

#endif
