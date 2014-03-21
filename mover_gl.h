#ifndef MOVER_GL_H
#define MOVER_GL_H

#ifdef GLFW
    #include "glfw_utils.h"
#else
    #include "GLES2/gl2.h"
    #include "egl_utils.h"
#endif

typedef struct mover_t
{
    char mover_type;

    // Program handle
    GLuint sphere_program;
    GLuint rectangle_program;

    // Locations
    GLint sphere_position_location;
    GLint sphere_center_location;
    GLint sphere_tex_coord_location;
    GLint sphere_color_location;
    GLint sphere_radius_location;

    // Locations
    GLint rectangle_position_location;
    GLint rectangle_center_location;
    GLint rectangle_tex_coord_location;
    GLint rectangle_color_location;

    // buffers
    GLuint vbo;
} mover_t;

void init_mover(mover_t *state);
void render_mover(float *center, float *gl_dims, float *color, mover_t *state);
void create_sphere_mover_program(mover_t *state);
void create_rectangle_mover_program(mover_t *state);
void draw_circle_mover(mover_t *state, float *center, float radius, float *color);
void draw_rectangle_mover(mover_t *state, float *center, float *color);
void create_mover_buffers(mover_t *state);

#endif
