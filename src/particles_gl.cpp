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

#include "particles_gl.h"

#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>

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

extern "C" void init_particles(particles_t *state, int screen_width, int screen_height)
{
    state->screen_width = screen_width;
    state->screen_height = screen_height;


    // Create circle buffers
    create_particle_buffers(state);

    // Create and set particle shaders
    // Also links particle program
    create_particle_shaders(state);
}

// Update coordinate of fluid points
extern "C" void render_particles(float *points, float diameter_pixels, int num_points, particles_t *state)
{
    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Orphan current buffer
    glBufferData(GL_ARRAY_BUFFER, 6*num_points*sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 6*num_points*sizeof(GLfloat), points, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    draw_particles(state, diameter_pixels, num_points);
}

void create_particle_buffers(particles_t *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

void create_particle_shaders(particles_t *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    compile_shader(vertexShader, "shaders/particle.vert");

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    compile_shader(fragmentShader, "shaders/particle.frag");

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
    state->radius_world_location = glGetUniformLocation(state->program, "radius_world");
    // Get pixel diameter location
    state->diameter_pixels_location = glGetUniformLocation(state->program, "diameter_pixels");
    // Get world to camera view matrix location
    state->view_matrix_location = glGetUniformLocation(state->program, "view_matrix_location");
    // Get camera to clip  projection matrix location
    state->proj_matrix_location = glGetUniformLocation(state->program, "proj_matrix_location");

    // Enable point size to be specified in the shader
    glEnable(GL_PROGRAM_POINT_SIZE);

//   GLfloat fSizes[2];
//   glGetFloatv(GL_POINT_SIZE_RANGE,fSizes);
//   printf("min: %f, max: %f\n", fSizes[0], fSizes[1]);
}

void draw_particles(particles_t *state, float diameter_pixels, int num_points)
{
    // Bind circle shader program
    glUseProgram(state->program);

    // Set radius uniform
    glUniform1f(state->radius_world_location, (GLfloat)diameter_pixels/state->screen_width);

    // Set pixel diameter uniform
    glUniform1f(state->diameter_pixels_location, (GLfloat)diameter_pixels);

    // Set view matrix
    glm::mat4 view = glm::lookAt(
        glm::vec3(1.2f, 1.2f, 1.2f),
        glm::vec3(0.0f, 0.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 1.0f)
    );
    glUniformMatrix4fv(state->view_matrix_location, 1, GL_FALSE, glm::value_ptr(view));

    // Set projection matrix
    float ratio = state->screen_width/state->screen_height;
    glm::mat4 proj = glm::perspective(45.0f, ratio, 0.0f, 10.0f);
    glUniformMatrix4fv(state->proj_matrix_location, 1, GL_FALSE, glm::value_ptr(proj));

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->color_location, 3, GL_FLOAT, GL_FALSE, 6*sizeof(GL_FLOAT),(void*)(3*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->color_location);

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw
    glDrawArrays(GL_POINTS, 0, num_points);
}
