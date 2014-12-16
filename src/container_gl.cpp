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
#include <iostream>
#include <stdlib.h>
#include "fluid.h"
#include "renderer.h"
#include "container_gl.h"

#include "ogl_utils.h"
#include "glfw_utils.h"
#include "world_gl.h"

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

void create_container_program(container_t *state)
{
    // Compile vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    compile_shader(vertex_shader, "shaders/container.vert");

    // Compile fragment shader
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    compile_shader(frag_shader, "shaders/container.frag");

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertex_shader);
    glAttachShader(state->program, frag_shader);

    // Link  program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Enable program
    glUseProgram(state->program);

    // Get position attribute location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get normal attribute location
    state->normal_location = glGetAttribLocation(state->program, "normal");
    // Get tex_coord attribute location
    state->tex_coord_location = glGetAttribLocation(state->program, "tex_coord");
    // Get color uniform location
    state->color_location = glGetUniformLocation(state->program, "color");
    // Get global matrix index
    state->global_matrix_index = glGetUniformBlockIndex(state->program, "GlobalMatrices");
    // Get global light index
    state->global_light_index = glGetUniformBlockIndex(state->program, "GlobalLight");

    // Setup buffers
    glBindVertexArray(state->vao);    
    size_t vert_size = 8*sizeof(GL_FLOAT);
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    glVertexAttribPointer(state->position_location, 3, GL_FLOAT, GL_FALSE, vert_size, 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->normal_location, 3, GL_FLOAT, GL_FALSE, vert_size, (void*)(3*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->normal_location);
    glVertexAttribPointer(state->tex_coord_location, 2, GL_FLOAT, GL_FALSE, vert_size, (void*)(6*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->tex_coord_location);

    glBindVertexArray(0);

    glUseProgram(0);
}

void create_container_buffers(container_t *state)
{
    // Generate VAO
    glGenVertexArrays(1, &state->vao);

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

void create_container_vertices(container_t *state)
{
    // vert x, y, z, normal x, y, z, tex_coord x, y
    float vertices[] = {
        // Floor
       -1.0, -1.0,  1.0, 0.0, 1.0, 0.0, 0.0, 0.0,
       -1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 1.0,
        1.0, -1.0,  1.0, 0.0, 1.0, 0.0, 1.0, 0.0,

       -1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 1.0,
        1.0, -1.0,  1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
        1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 1.0, 1.0,
        //right wall
        1.0, -1.0,  1.0, -1.0, 0.0, 0.0, 1.0, 0.0,
        1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 1.0,  1.0,  -1.0, 0.0, 0.0, 1.0, 1.0, 

        1.0, -1.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0,
        1.0, 1.0,  1.0,  -1.0, 0.0, 0.0, 1.0, 1.0,
        1.0, 1.0, -1.0,  -1.0, 0.0, 0.0, 0.0, 1.0,
        // Back
       -1.0, 1.0,  -1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
        1.0, 1.0,  -1.0, 0.0, 0.0, 1.0, 1.0, 1.0,
        1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0,

        1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0,
       -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
       -1.0,  1.0, -1.0, 0.0, 0.0, 1.0, 0.0, 1.0,
        // Left wall
       -1.0, -1.0,  1.0, 1.0, 0.0, 0.0, 1.0, 0.0,
       -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
       -1.0, 1.0,  1.0,  1.0, 0.0, 0.0, 1.0, 1.0,

       -1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
       -1.0, 1.0,  1.0,  1.0, 0.0, 0.0, 1.0, 1.0,
       -1.0, 1.0, -1.0,  1.0, 0.0, 0.0, 0.0, 1.0,
    };

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 8*24*sizeof(GLfloat), vertices, GL_STATIC_DRAW);
}

void init_container(container_t *state, int screen_width, int screen_height)
{
    // Set screen with/height in pixels
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Generate buffers 
    create_container_buffers(state);

    // Create program
    create_container_program(state);

    // Set verticies
    create_container_vertices(state);
}

void render_container(container_t *state)
{
    // Setup program
    glUseProgram(state->program);

    // set color uniform
    float color[] = {1.0, 1.0, 1.0, 1.0};
    glUniform4fv(state->color_location, 1, color);

    // Set uniform binding
    glUniformBlockBinding(state->program, state->global_matrix_index, g_GlobalMatricesBindingIndex);
    glUniformBlockBinding(state->program, state->global_light_index, g_GlobalLightBindingIndex);

    // Bind VAO
    glBindVertexArray(state->vao);   

    // Draw container
    glDrawArrays(GL_TRIANGLES, 0, 24);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

}
