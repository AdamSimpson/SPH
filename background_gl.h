#ifndef BACKGROUND_GL_H
#define BACKGROUND_GL_H

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

typedef struct background_t {
    GLuint program;

    // Program locations
    GLint position_location;
    GLint tex_coord_location;
    GLint tex_location;

    // Uniforms
    GLuint tex_uniform;

    // Vertex buffer
    GLuint vbo;

    // Element buffer
    GLuint ebo;

} background_t;

void create_backround_program(background_t *state);
void create_background_buffers(background_t *state);
void create_background_vertices(background_t *state);
void create_background_texture(background_t *state);
void init_background(background_t *state);
void draw_background(background_t *state);

#endif
