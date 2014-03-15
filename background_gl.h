#ifndef BACKGROUND_GL_H
#define BACKGROUND_GL_H

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

typedef struct {
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

} BACKGROUND_T;

void create_backround_program(BACKGROUND_T *state);
void create_background_buffers(BACKGROUND_T *state);
void create_background_vertices(BACKGROUND_T *state);
void create_background_texture(BACKGROUND_T *state);
void init_background(BACKGROUND_T *state);
void draw_background(BACKGROUND_T *state);

#endif
