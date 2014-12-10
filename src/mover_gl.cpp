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

#ifdef __cplusplus
extern "C" {
#endif
#include "ogl_utils.h"
#include "glfw_utils.h"
#include "world_gl.h"
#ifdef __cplusplus
}
#endif

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

void init_mover(mover_t *state, int screen_width, int screen_height)
{
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Create VAO
    create_mover_buffers(state);

    // Create shader programs
    create_sphere_mover_program(state);
}

void create_mover_buffers(mover_t *state)
{
    // VAO
    glGenVertexArrays(1, &state->vao);
}


// Update coordinates of point mover and then render
void render_mover(float *center, float radius, float *color, mover_t *state)
{
    draw_circle_mover(state, center, radius, color);
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

    // Compile geometry shader
    GLuint geometryShader = glCreateShader(GL_GEOMETRY_SHADER);
    compile_shader(geometryShader, "shaders/mover_circle.geom");

    // Create shader program
    state->sphere_program = glCreateProgram();
    glAttachShader(state->sphere_program, vertexShader);
    glAttachShader(state->sphere_program, fragmentShader); 
    glAttachShader(state->sphere_program, geometryShader);

    // Link and use program
    glLinkProgram(state->sphere_program);
    show_program_log(state->sphere_program);

    // Get color location
    state->sphere_color_location = glGetUniformLocation(state->sphere_program, "color");
    // Get radius location
    state->sphere_radius_location = glGetUniformLocation(state->sphere_program, "sphereRadius");
    // Get center location
    state->sphere_center_location = glGetUniformLocation(state->sphere_program, "center");
    // Get global matrix index
    state->global_matrix_index = glGetUniformBlockIndex(state->sphere_program, "GlobalMatrices");
    // Get global light index
    state->global_light_index = glGetUniformBlockIndex(state->sphere_program, "GlobalLight");
}

void draw_circle_mover(mover_t *state, float *center, float radius, float *color)
{
    // Bind sphere shader program
    glUseProgram(state->sphere_program);

    // Enable VAO
    glBindVertexArray(state->vao);

    // set radius uniform
    glUniform1f(state->sphere_radius_location, radius);

    // set color uniform
    glUniform3fv(state->sphere_color_location, 1, color);

    // set center uniform
    glUniform3fv(state->sphere_center_location, 1, center);

    // Set uniform binding
    glUniformBlockBinding(state->sphere_program, state->global_matrix_index, g_GlobalMatricesBindingIndex);
    glUniformBlockBinding(state->sphere_program, state->global_light_index, g_GlobalLightBindingIndex);

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw
    glDrawArrays(GL_POINTS, 0, 1);
  
    glDisable(GL_BLEND);
}
