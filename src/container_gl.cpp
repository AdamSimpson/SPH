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
#include <stdlib.h>
#include "fluid.h"
#include "renderer.h"
#include "container_gl.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "ogl_utils.h"
#include "glfw_utils.h"
#ifdef __cplusplus
}
#endif

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

    // Get position attribute location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get color uniform location
    state->color_location = glGetUniformLocation(state->program, "color");
    // Get world to camera view matrix location
    state->view_matrix_location = glGetUniformLocation(state->program, "view");
    // Get camera to clip  projection matrix location
    state->proj_matrix_location = glGetUniformLocation(state->program, "proj");
}

void create_container_buffers(container_t *state)
{
    // VAO is required for OpenGL 3+ when using VBO I believe
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);

    // Generate element buffer
    glGenBuffers(1, &state->ebo);
}

void create_container_vertices(container_t *state)
{
    float vertices[] = {
         // Full screen vertices
        0.0, 0.0, -1.0,
        0.0, 0.0, 1.0,
        1.0, 0.0, -1.0,
        1.0, 0.0, 1.0
    };

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 3*4*3*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

    // Elements
    GLubyte elements[] = {
        1, 0, 2,
        3, 1, 2
    };

    // Set buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state->ebo);
    // Fill buffer
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2*3*sizeof(GLubyte), elements, GL_STATIC_DRAW);
}

void init_container(container_t *state, int screen_width, int screen_height)
{
    // Set screen with/height in pixels
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Create program
    create_container_program(state);

    // Generate buffers 
    create_container_buffers(state);  

    // Set verticies
    create_container_vertices(state);
}

void draw_container(container_t *state)
{
    // Setup program
    glUseProgram(state->program);

    // Setup buffers
    size_t vert_size = 3*sizeof(GL_FLOAT);
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, vert_size, 0);
    glEnableVertexAttribArray(state->position_location);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state->ebo);

    // set color uniform
    float color[] = {1.0, 1.0, 1.0, 1.0};
    glUniform4fv(state->color_location, 1, color);

   // Set view matrix
    glm::mat4 view = glm::lookAt(
        glm::vec3(0.0f, 0.2f, 0.2f), // Eye position
        glm::vec3(0.0f, 0.0f, 0.0f), // Looking at
        glm::vec3(0.0f, 1.0f, 0.0f)  // Up
    );

    glUniformMatrix4fv(state->view_matrix_location, 1, GL_FALSE, glm::value_ptr(view));

    // Set projection matrix
    float ratio = (float)state->screen_width/(float)state->screen_height;
    glm::mat4 proj = glm::perspective(45.0f, ratio, 1.0f, 10.0f);
    glUniformMatrix4fv(state->proj_matrix_location, 1, GL_FALSE, glm::value_ptr(proj));

    // Disable Blend
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw container
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);
}
