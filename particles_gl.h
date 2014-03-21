#ifndef PARTICLES_GL_H
#define PARTICLES_GL_H

#ifdef GLFW
    #include "glfw_utils.h"
#else
    #include "GLES2/gl2.h"
    #include "egl_utils.h"
#endif

typedef struct particles_t {
    // Program handle
    GLuint program;

    // Locations
    GLint position_location;
    GLint color_location;
    GLint diameter_pixels_location;
    GLint radius_world_location;

    // Screen dimensions
    int screen_width;
    int screen_height;

    // buffers
    GLuint vbo;
} particles_t;

void init_particles(particles_t *state, int screen_width, int screen_height);
void render_particles(float *points, float diameter_pixels, int num_points, particles_t *state);
void create_particle_shaders(particles_t *state);
void draw_particles(particles_t *state, float diameter_pixels, int num_points);
void create_particle_buffers(particles_t *state);

#endif
