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

void init_world(world_t *state, int screen_width, int screen_height)
{
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Set global matrices
    glGenBuffers(1, &g_GlobalMatricesUBO);
    glBindBuffer(GL_UNIFORM_BUFFER, g_GlobalMatricesUBO);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::mat4) * 2, NULL, GL_STREAM_DRAW);
    // view matric
    glm::mat4 view = glm::lookAt(
        glm::vec3(0.0f, 0.2f, 0.2f), // Eye position
        glm::vec3(0.0f, 0.0f, 0.0f), // Looking at
        glm::vec3(0.0f, 1.0f, 0.0f)  // Up
    );
    // Set projection matrix
    float ratio = (float)state->screen_width/(float)state->screen_height;
    glm::mat4 proj = glm::perspective(45.0f, ratio, 1.0f, 10.0f);
    glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(glm::mat4), glm::value_ptr(view));
    glBufferSubData(GL_UNIFORM_BUFFER, sizeof(glm::mat4), sizeof(glm::mat4), glm::value_ptr(proj));

    glBindBufferRange(GL_UNIFORM_BUFFER, g_GlobalMatricesBindingIndex, g_GlobalMatricesUBO, 0, sizeof(glm::mat4) * 2);
}
