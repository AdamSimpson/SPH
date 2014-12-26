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
#include "camera.hpp"

#include "ogl_utils.h"
#include "gl.hpp"

#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include "glm/gtx/norm.hpp"

GLuint Camera::g_GlobalMatricesUBO;
GLuint Camera::g_GlobalLightUBO;

Camera::Camera(GL& gl): gl(gl)
{
    // Camera position rotations are based upon
    this->position_default = glm::vec3(0.0f, 0.0f, 0.8f);

    // current camera position
    this->position = glm::vec3(0.0f, 0.0f, 0.8f);

    // Camera look at point
    this->look_at = glm::vec3(0.0, -0.5f, -0.7f);

    // perspective zoom factor
    this->zoom_factor = 1.0f;

    // Set global matrices
    glGenBuffers(1, &Camera::g_GlobalMatricesUBO);

    // Create and Allocate buffer storage
    glBindBuffer(GL_UNIFORM_BUFFER, Camera::g_GlobalMatricesUBO);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::mat4) * 2, NULL, GL_STREAM_DRAW);

    // Attach to binding index
    glBindBufferRange(GL_UNIFORM_BUFFER, Camera::g_GlobalMatricesBindingIndex, Camera::g_GlobalMatricesUBO, 0, sizeof(glm::mat4) * 2);

    glGenBuffers(1, &Camera::g_GlobalLightUBO);

    // Create and Allocate buffer storage
    glBindBuffer(GL_UNIFORM_BUFFER, Camera::g_GlobalLightUBO);
    glBufferData(GL_UNIFORM_BUFFER, 4*sizeof(glm::vec4) + sizeof(float), NULL, GL_STREAM_DRAW);

    // Attach to binding index
    glBindBufferRange(GL_UNIFORM_BUFFER, Camera::g_GlobalLightBindingIndex, Camera::g_GlobalLightUBO, 0, 4*sizeof(glm::vec4) + sizeof(float));

    // Set all uniform values
    this->update_view();

    // Unbind
    glBindBuffer(GL_UNIFORM_BUFFER, 0);

}

// Rotates the default camera  position degrees around the y axis
void Camera::rotate_yaw(float degrees)
{
    // We wish to rotate the vector from the look at point to the camera position
    glm::vec3 look_to_camera_default = this->position_default - this->look_at;

    // Rotate vector
    float rads = degrees * M_PI / 180.0f;
    glm::vec3 look_to_camera_new = glm::rotateY(look_to_camera_default, rads);

    // Regain camera vector by adding in look_at position
    this->position = look_to_camera_new + look_at;
}

// Rotates the default camera position degrees around the x axis
void Camera::rotate_pitch(float degrees)
{
    // We wish to rotate the vector from the look at point to the camera position
    glm::vec3 look_to_camera_default = this->position_default - this->look_at;

    // Rotate vector
    float rads = degrees * M_PI / 180.0f;
    glm::vec3 look_to_camera_new = glm::rotateX(look_to_camera_default, rads);

    // Regain camera vector by adding in look_at position
    this->position = look_to_camera_new + this->look_at;
}

void Camera::rotate_yaw_pitch(float degrees_yaw, float degrees_pitch)
{
    // We wish to rotate the vector from the look at point to the camera position
    glm::vec3 look_to_camera_default = this->position_default - this->look_at;

    // Rotate vector
    float rad_pitch = degrees_pitch * M_PI / 180.0f;
    float rad_yaw = degrees_yaw * M_PI / 180.0f;
    glm::vec3 look_to_camera_new = glm::rotateY(look_to_camera_default, rad_yaw);
    look_to_camera_new = glm::rotateX(look_to_camera_new, rad_pitch);

    // Regain camera vector by adding in look_at position
    this->position = look_to_camera_new + look_at;
}

// Move along view vector
void Camera::move_along_view(float dx)
{
    // Get vector from look_at position to camera
    glm::vec3 look_to_camera_old = this->position - this->look_at;
    
    // Get normal along look_at to camera vector
    float norm = glm::l1Norm(look_to_camera_old);
    glm::vec3 look_to_camera_norm = look_to_camera_old/norm;
    
    // subtract dx along the vector
    glm::vec3 look_to_camera_new = look_to_camera_old - look_to_camera_norm*dx;

    // Regain camera vector by adding in look_at position
    this->position = look_to_camera_new + this->look_at;

    // We have to update the default as well to account for move

    // Get vector from look_at position to camera
    look_to_camera_old = this->position_default - this->look_at;

    // Get normal along look_at to camera vector
    norm = glm::l1Norm(look_to_camera_old);
    look_to_camera_norm = look_to_camera_old/norm;

    // subtract dx along the vector
    look_to_camera_new = look_to_camera_old - look_to_camera_norm*dx;

    // Regain camera vector by adding in look_at position
    this->position_default = look_to_camera_new + this->look_at;
}

// Zoom persective...Not a huge fan of this method
// But it's easier to deal with
void Camera::zoom_view(float d_zoom_factor)
{
    this->zoom_factor += d_zoom_factor;
}

void Camera::update_view()
{
    // Bind UBO
    glBindBuffer(GL_UNIFORM_BUFFER, Camera::g_GlobalMatricesUBO);
    
    // Update view matrix
    glm::mat4 view = glm::lookAt(
        this->position,
        this->look_at, // Looking at
        glm::vec3(0.0f, 1.0f, 0.0f)  // Up
    );
    // Set projection matrix
    float ratio = (float)this->gl.get_screen_width()/(float)this->gl.get_screen_height();
    // 1.22 radians ~ 70 degrees
    glm::mat4 proj = glm::perspective(this->zoom_factor*1.22f, ratio, 0.7f, 10.0f);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4), glm::value_ptr(view));
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4), glm::value_ptr(proj));

    // Create and Allocate buffer storage
    glBindBuffer(GL_UNIFORM_BUFFER, g_GlobalLightUBO);

    // Update light attributes
    glm::vec4 worldSpacePos = glm::vec4(0.3, 0.5, -0.4, 1.0);
    glm::vec4 cameraSpacePos     = view*worldSpacePos;
    glm::vec4 intensity     = glm::vec4(0.8, 0.8, 0.8, 1.0);
    glm::vec4 ambientIntensity   = glm::vec4(0.1, 0.1, 0.1, 1.0);
    float attenuation       = 1.0f;

    // Buffer uniform data
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::vec4), glm::value_ptr(worldSpacePos));
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::vec4), sizeof(glm::vec4), glm::value_ptr(cameraSpacePos));
    glBufferSubData(GL_UNIFORM_BUFFER, 2*sizeof(glm::vec4), sizeof(glm::vec4), glm::value_ptr(intensity));
    glBufferSubData(GL_UNIFORM_BUFFER, 3*sizeof(glm::vec4), sizeof(glm::vec4), glm::value_ptr(ambientIntensity));
    glBufferSubData(GL_UNIFORM_BUFFER, 4*sizeof(glm::vec4), sizeof(float), &attenuation);

    // Unbind
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
}
