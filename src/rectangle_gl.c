/*
The MIT License (MIT)

Copyright (c) 2014 Adam Simpson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "rectangle_gl.h"
#include "ogl_utils.h"
#include "glfw_utils.h"

void init_rectangle(rectangle_t *state)
{
    // Create rectangle buffers
    create_rectangle_buffers(state);

    // Create shader programs
    create_rectangle_program(state);
}

// Update coordinates of point rectangle and then draw
void render_rectangle(rectangle_t *state, float *center, float *gl_dims, float *color)
{
     // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    float center_x = center[0];
    float center_y = center[1];

    float half_width = gl_dims[0]*0.5;
    float half_height = gl_dims[1]*0.5;
    // Triangle strip Verticies for rectangle that contains circle
    float verts[2*4];

    // Verticies = gl_x, gl_y, tex_x, tex_y

    // Bottom left
    verts[0] = center_x - half_width; 
    verts[1] = center_y - half_height;
    // Top left
    verts[2] = center_x - half_width;
    verts[3] = center_y + half_height;
    // Bottom right
    verts[4] = center_x + half_width;
    verts[5] = center_y - half_height;
    // Top right
    verts[6] = center_x + half_width;
    verts[7] = center_y + half_height;

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 2*4*sizeof(GLfloat), verts, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Draw rectangle
    draw_rectangle(state, center, color);
}

void create_rectangle_buffers(rectangle_t *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    #ifndef RASPI
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

// Compile rectnagle program
void create_rectangle_program(rectangle_t *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
      compile_shader(vertexShader, "SPH/shaders/rectangle_es.vert");
    #else
      compile_shader(vertexShader, "shaders/rectangle.vert");
    #endif

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
      compile_shader(fragmentShader, "SPH/shaders/rectangle_es.frag");
    #else
      compile_shader(fragmentShader, "shaders/rectangle.frag");
    #endif

    // Create shader program
    state->rectangle_program = glCreateProgram();
    glAttachShader(state->rectangle_program, vertexShader);
    glAttachShader(state->rectangle_program, fragmentShader);

    // Link and use program
    glLinkProgram(state->rectangle_program);
    show_program_log(state->rectangle_program);

    // Get rectangle location
    state->color_location = glGetUniformLocation(state->rectangle_program, "color");
    state->position_location = glGetAttribLocation(state->rectangle_program, "position");
}

void draw_rectangle(rectangle_t *state, float *center, float *color)
{
    // Bind rectangle shader program
    glUseProgram(state->rectangle_program);

    // set color uniform
    glUniform4fv(state->color_location, 1, color);

    // Bind buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Setup verticies
    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

