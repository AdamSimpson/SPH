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

    // Create shader programs
    create_sphere_mover_program(state);
    create_rectangle_mover_program(state);
}

// Update coordinates of point mover and then render
void render_mover(float *center, float *gl_dims, float *color, MOVER_T *state)
{
     // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    float center_x = center[0];
    float center_y = center[1];

    float half_width = gl_dims[0]*0.5;
    float half_height = gl_dims[1]*0.5;
    // Triangle strip Verticies for rectangle that contains circle
    float verts[4*4];

    // Verticies = gl_x, gl_y, tex_x, tex_y

    // Bottom left
    verts[0] = center_x - half_width; 
    verts[1] = center_y - half_height;
    verts[2] = -1.0;
    verts[3] = -1.0;
    // Top left
    verts[4] = center_x - half_width;
    verts[5] = center_y + half_height;
    verts[6] = -1.0;
    verts[7] = 1.0;
    // Bottom right
    verts[8] = center_x + half_width;
    verts[9] = center_y - half_height;
    verts[10] = 1.0;
    verts[11] = -1.0;
    // Top right
    verts[12] = center_x + half_width;
    verts[13] = center_y + half_height;
    verts[14] = 1.0;
    verts[15] = 1.0;

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 4*4*sizeof(GLfloat), verts, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    if(state->mover_type == SPHERE_MOVER) {
        float radius = half_width;
        draw_circle_mover(state, center, radius, color);
    }
    else if(state->mover_type == RECTANGLE_MOVER)
        draw_rectangle_mover(state, center, color);
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

// Compile sphere program
void create_sphere_mover_program(MOVER_T *state)
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
    state->sphere_program = glCreateProgram();
    glAttachShader(state->sphere_program, vertexShader);
    glAttachShader(state->sphere_program, fragmentShader); 

    // Link and use program
    glLinkProgram(state->sphere_program);
    show_program_log(state->sphere_program);

    // Get position location
    state->sphere_position_location = glGetAttribLocation(state->sphere_program, "position");
    state->sphere_tex_coord_location = glGetAttribLocation(state->sphere_program, "tex_coord");

    // Get color location
    state->sphere_color_location = glGetUniformLocation(state->sphere_program, "color");
    // Get radius location
    state->sphere_radius_location = glGetUniformLocation(state->sphere_program, "radius");
    // Get center location
    state->sphere_center_location = glGetUniformLocation(state->sphere_program, "center");

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

// Compile rectnagle program
void create_rectangle_mover_program(MOVER_T *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef GLES
      compile_shader(vertexShader, "SPH/mover_rectangle_es.vert");
    #else
      compile_shader(vertexShader, "mover_rectangle.vert");
    #endif

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef GLES
      compile_shader(fragmentShader, "SPH/mover_rectangle_es.frag");
    #else
      compile_shader(fragmentShader, "mover_rectangle.frag");
    #endif

    // Create shader program
    state->rectangle_program = glCreateProgram();
    glAttachShader(state->rectangle_program, vertexShader);
    glAttachShader(state->rectangle_program, fragmentShader);

    // Link and use program
    glLinkProgram(state->rectangle_program);
    show_program_log(state->rectangle_program);

    // Get position location
    state->rectangle_position_location = glGetAttribLocation(state->rectangle_program, "position");
    state->rectangle_tex_coord_location = glGetAttribLocation(state->rectangle_program, "tex_coord");

    // Get rectangle location
    state->rectangle_color_location = glGetUniformLocation(state->rectangle_program, "color");
    // Get center location
    state->rectangle_center_location = glGetUniformLocation(state->rectangle_program, "center");
}

void draw_circle_mover(MOVER_T *state, float *center, float radius, float *color)
{

    // Bind sphere shader program
    glUseProgram(state->sphere_program);

    // set radius uniform
    glUniform1f(state->sphere_radius_location, radius);

    // set color uniform
    glUniform3fv(state->sphere_color_location, 1, color);

    // set center uniform
    glUniform2fv(state->sphere_center_location, 1, center);

    // Bind buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Setup verticies
    glVertexAttribPointer(state->sphere_position_location, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GL_FLOAT), 0);
    glVertexAttribPointer(state->sphere_tex_coord_location, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GL_FLOAT), (void*)(2*sizeof(GLfloat)));
    glEnableVertexAttribArray(state->sphere_position_location);

    // Draw
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void draw_rectangle_mover(MOVER_T *state, float *center, float *color)
{

    // Bind rectangle shader program
    glUseProgram(state->rectangle_program);

    // set color uniform
    glUniform3fv(state->rectangle_color_location, 1, color);

    // set center uniform
    glUniform2fv(state->rectangle_center_location, 1, center);

    // Bind buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Setup verticies
    glVertexAttribPointer(state->rectangle_position_location, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GL_FLOAT), 0);
    glVertexAttribPointer(state->rectangle_position_location, 2, GL_FLOAT, GL_FALSE, 4*sizeof(GL_FLOAT), (void*)(2*sizeof(GLfloat)));
    glEnableVertexAttribArray(state->rectangle_position_location);

    // Draw
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

