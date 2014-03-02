#ifndef PARTICLES_GL_H
#define PARTICLES_GL_H

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
} PARTICLES_T;

void init_particles(PARTICLES_T *state);
void update_particles(float *points, float radius, int num_points, PARTICLES_T *state);
void create_particle_shaders(PARTICLES_T *state);
void draw_particles(PARTICLES_T *state, float radius, int num_points);
void create_particle_buffers(PARTICLES_T *state);

#endif
