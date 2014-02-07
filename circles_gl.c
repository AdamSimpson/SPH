#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "circles_gl.h"
#include "ogl_utils.h"

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

void init_circles(CIRCLE_T *state)
{
    // Create circle buffers
    create_buffers(state);

    // Create and set circle shaders
    // Also links circle program
    create_shaders(state);
}

// Update coordinates of point mover
void update_mover_point(float *point, float radius, CIRCLE_T *state)
{
     // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 5*sizeof(GLfloat), point, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    draw_circle_mover(state, radius);
}

// Update coordinate of fluid points
void update_points(float *points, float radius, int num_points, CIRCLE_T *state)
{
    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 5*num_points*sizeof(GLfloat), points, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    draw_circles(state, radius, num_points);
}

void create_buffers(CIRCLE_T *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    #ifndef GLES
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

void create_shaders(CIRCLE_T *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef GLES
      compile_shader(vertexShader, "SPH/particle_es.vert");
    #else
      compile_shader(vertexShader, "particle.vert");
    #endif

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef GLES
      compile_shader(fragmentShader, "SPH/particle_es.frag");
    #else
      compile_shader(fragmentShader, "particle.frag");
    #endif

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertexShader);
    glAttachShader(state->program, fragmentShader); 

    // Link and use program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Get position location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get tex_coord location
    state->color_location = glGetAttribLocation(state->program, "color");
    // Get radius location
    state->radius_location = glGetUniformLocation(state->program, "radius");

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Enable point size to be specified in the shader
    #ifndef GLES
    glEnable(GL_PROGRAM_POINT_SIZE);
    #endif

}

void draw_circle_mover(CIRCLE_T *state, float radius)
{
    // Bind circle shader program
    glUseProgram(state->program);

    // set radius uniform
    glUniform1f(state->radius_location, (GLfloat)radius);

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->color_location, 3, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT),(void*)(2*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->color_location);

    // Draw
    glDrawArrays(GL_POINTS, 0, 1);

}

void draw_circles(CIRCLE_T *state, float radius, int num_points)
{
    // Bind circle shader program
    glUseProgram(state->program);

    // set radius uniform
    glUniform1f(state->radius_location, (GLfloat)radius);

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->color_location, 3, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT),(void*)(2*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->color_location);

    // Draw
    glDrawArrays(GL_POINTS, 0, num_points);
}
