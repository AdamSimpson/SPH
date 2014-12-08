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
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "world_gl.h"

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
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtx/norm.hpp"

void init_world(world_t *state, int screen_width, int screen_height)
{
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    state->max_degrees_rotate = 40.0f;

    // eye position rotations are based upon
    state->eye_position_default[0] = 0.0f;
    state->eye_position_default[1] = 0.0f;
    state->eye_position_default[2] = 0.8f;

    // current eye position
    state->eye_position[0] = 0.0f;
    state->eye_position[1] = 0.0f;
    state->eye_position[2] = 0.8f;

    state->look_at[0] = 0.0f;
    state->look_at[1] = -0.5;
    state->look_at[2] = -0.7;

    state->zoom_factor = 1.0f;

    // Set global matrices
    glGenBuffers(1, &g_GlobalMatricesUBO);

    // Create and Allocate buffer storage
    glBindBuffer(GL_UNIFORM_BUFFER, g_GlobalMatricesUBO);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::mat4) * 2, NULL, GL_STREAM_DRAW);

    // Update buffers
    update_view(state);

    // Attach to binding index
    glBindBuffer(GL_UNIFORM_BUFFER, g_GlobalMatricesUBO);
    glBindBufferRange(GL_UNIFORM_BUFFER, g_GlobalMatricesBindingIndex, g_GlobalMatricesUBO, 0, sizeof(glm::mat4) * 2);
}

// Rotates the default eye position degrees around the y axis
void rotate_camera_yaw(world_t *state, float degrees)
{
    glm::vec3 eye_default = glm::vec3(state->eye_position_default[0], state->eye_position_default[1], state->eye_position_default[2]);
    glm::vec3 look_at = glm::vec3(state->look_at[0], state->look_at[1], state->look_at[2]);

    // We wish to rotate the vector from the look at point to the eye position
    glm::vec3 look_to_eye_default = eye_default - look_at;

    // Rotate vector
    float rads = degrees * M_PI / 180.0f;
    glm::vec3 look_to_eye_new = glm::rotateY(look_to_eye_default, rads);

    // Regain eye vector by adding in look_at position
    glm::vec3 eye_new = look_to_eye_new + look_at;

    state->eye_position[0] = eye_new[0];
    state->eye_position[1] = eye_new[1];
    state->eye_position[2] = eye_new[2];
}

// Rotates the default eye position degrees around the x axis
void rotate_camera_pitch(world_t *state, float degrees)
{
    glm::vec3 eye_default = glm::vec3(state->eye_position_default[0], state->eye_position_default[1], state->eye_position_default[2]);
    glm::vec3 look_at = glm::vec3(state->look_at[0], state->look_at[1], state->look_at[2]);

    // We wish to rotate the vector from the look at point to the eye position
    glm::vec3 look_to_eye_default = eye_default - look_at;

    // Rotate vector
    float rads = degrees * M_PI / 180.0f;
    glm::vec3 look_to_eye_new = glm::rotateX(look_to_eye_default, rads);

    // Regain eye vector by adding in look_at position
    glm::vec3 eye_new = look_to_eye_new + look_at;

    state->eye_position[0] = eye_new[0];
    state->eye_position[1] = eye_new[1];
    state->eye_position[2] = eye_new[2];
}

void rotate_camera_yaw_pitch(world_t *state, float degrees_yaw, float degrees_pitch)
{
    glm::vec3 eye_default = glm::vec3(state->eye_position_default[0], state->eye_position_default[1], state->eye_position_default[2]);
    glm::vec3 look_at = glm::vec3(state->look_at[0], state->look_at[1], state->look_at[2]);

    // We wish to rotate the vector from the look at point to the eye position
    glm::vec3 look_to_eye_default = eye_default - look_at;

    // Rotate vector
    float rad_pitch = degrees_pitch * M_PI / 180.0f;
    float rad_yaw = degrees_yaw * M_PI / 180.0f;
    glm::vec3 look_to_eye_new = glm::rotateY(look_to_eye_default, rad_yaw);
    look_to_eye_new = glm::rotateX(look_to_eye_new, rad_pitch);

    // Regain eye vector by adding in look_at position
    glm::vec3 eye_new = look_to_eye_new + look_at;

    state->eye_position[0] = eye_new[0];
    state->eye_position[1] = eye_new[1];
    state->eye_position[2] = eye_new[2];
}

// Move along view vector
void move_along_view(world_t *state, float dx)
{
    glm::vec3 eye_old = glm::vec3(state->eye_position[0], state->eye_position[1], state->eye_position[2]);
    glm::vec3 look_at = glm::vec3(state->look_at[0], state->look_at[1], state->look_at[2]);
    
    // Get vector from look_at position to eye
    glm::vec3 look_to_eye_old = eye_old - look_at;
    
    // Get normal along look_at to eye vector
    float norm = glm::l1Norm(look_to_eye_old);
    glm::vec3 look_to_eye_norm = look_to_eye_old/norm;
    
    // subtract dx along the vector
    glm::vec3 look_to_eye_new = look_to_eye_old - look_to_eye_norm*dx;

    // Regain eye vector by adding in look_at position
    glm::vec3 eye_new = look_to_eye_new + look_at;

    state->eye_position[0] = eye_new[0];
    state->eye_position[1] = eye_new[1];
    state->eye_position[2] = eye_new[2];

    // We have to update the default as well to account for move
    eye_old = glm::vec3(state->eye_position_default[0], state->eye_position_default[1], state->eye_position_default[2]);
    look_at = glm::vec3(state->look_at[0], state->look_at[1], state->look_at[2]);

    // Get vector from look_at position to eye
    look_to_eye_old = eye_old - look_at;

    // Get normal along look_at to eye vector
    norm = glm::l1Norm(look_to_eye_old);
    look_to_eye_norm = look_to_eye_old/norm;

    // subtract dx along the vector
    look_to_eye_new = look_to_eye_old - look_to_eye_norm*dx;

    // Regain eye vector by adding in look_at position
    eye_new = look_to_eye_new + look_at;

    state->eye_position_default[0] = eye_new[0];
    state->eye_position_default[1] = eye_new[1];
    state->eye_position_default[2] = eye_new[2];
}

// Zoom persective...Not a huge fan of this method
// But it's easier to deal with
void zoom_view(world_t *state, float d_zoom_factor)
{
    state->zoom_factor += d_zoom_factor;
}

void update_view(world_t *state)
{
    // Bind UBO
    glBindBuffer(GL_UNIFORM_BUFFER, g_GlobalMatricesUBO);
    
    // Update view matrix
    glm::mat4 view = glm::lookAt(
        glm::vec3(state->eye_position[0], state->eye_position[1], state->eye_position[2]), // Eye position
        glm::vec3(state->look_at[0], state->look_at[1], state->look_at[2]), // Looking at
        glm::vec3(0.0f, 1.0f, 0.0f)  // Up
    );
    // Set projection matrix
    float ratio = (float)state->screen_width/(float)state->screen_height;
    // 1.22 radians ~ 70 degrees
    glm::mat4 proj = glm::perspective(state->zoom_factor*1.22f, ratio, 0.7f, 10.0f);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4), glm::value_ptr(view));
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4), glm::value_ptr(proj));

    // Unbind
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}
