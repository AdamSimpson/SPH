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
#include "mover_gl.h"
#include "ogl_utils.h"

#include "glfw_utils.h"

void init_mover(mover_t *state)
{
    // Create circle buffers
    create_mover_buffers(state);

    // Create shader programs
    create_sphere_mover_program(state);
}

// Update coordinates of point mover and then render
void render_mover(float *center, float *gl_dims, float *color, mover_t *state)
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

    float radius = half_width;
    draw_circle_mover(state, center, radius, color);
}

void create_mover_buffers(mover_t *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

// Compile sphere program
void create_sphere_mover_program(mover_t *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    compile_shader(vertexShader, "shaders/mover_circle.vert");

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    compile_shader(fragmentShader, "shaders/mover_circle.frag");

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
}

void draw_circle_mover(mover_t *state, float *center, float radius, float *color)
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

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}
