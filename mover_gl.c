#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "mover_gl.h"
#include "ogl_utils.h"

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

void init_mover(MOVER_T *state)
{
    // Create circle buffers
    create_mover_buffers(state);

    // Create and set circle shaders
    // Also links circle program
    create_mover_shaders(state);
}

// Update coordinates of point mover
void update_mover(float *center, float *gl_dims, float *color, MOVER_T *state)
{
     // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    float center_x = center[0];
    float center_y = center[1];

    // Triangle strip Verticies for rectangle that contains circle
    float verts[8];

    // Bottom left
    verts[0] = center_x - gl_dims[0]; 
    verts[1] = center_y - gl_dims[1];
    // Top left
    verts[2] = center_x - gl_dims[0];
    verts[3] = center_y + gl_dims[1];
    // Bottom right
    verts[4] = center_x + gl_dims[0];
    verts[5] = center_y - gl_dims[1];
    // Top right
    verts[6] = center_x + gl_dims[0];
    verts[7] = center_y + gl_dims[1];

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 8*sizeof(GLfloat), verts, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    float radius = gl_dims[0];

    draw_circle_mover(state, center, radius, color);
}

void create_mover_buffers(MOVER_T *state)
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

void create_mover_shaders(MOVER_T *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef GLES
      compile_shader(vertexShader, "SPH/mover_circle_es.vert");
    #else
      compile_shader(vertexShader, "mover_circle.vert");
    #endif

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef GLES
      compile_shader(fragmentShader, "SPH/mover_circle_es.frag");
    #else
      compile_shader(fragmentShader, "mover_circle.frag");
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

    // Get tex_coorad location
    state->color_location = glGetUniformLocation(state->program, "color");
    // Get radius location
    state->radius_location = glGetUniformLocation(state->program, "radius");
    // Get center location
    state->center_location = glGetUniformLocation(state->program, "center");

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Enable point size to be specified in the shader
    #ifndef GLES
    glEnable(GL_PROGRAM_POINT_SIZE);
    #endif
}

void draw_circle_mover(MOVER_T *state, float *center, float radius, float *color)
{

    // Bind circle shader program
    glUseProgram(state->program);

    // set radius uniform
    glUniform1f(state->radius_location, radius);

    // set color uniform
    glUniform3fv(state->color_location, 1, color);

    // set center uniform
    glUniform2fv(state->center_location, 1, center);

    // Bind buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Setup verticies
    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);

    // Draw
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
