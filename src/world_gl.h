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

#ifndef WORLD_GL_H
#define WORLD_GL_H

#include "glfw_utils.h"
#include "structs.h"

// Binding index for global uniform matrices
static const int g_GlobalMatricesBindingIndex = 0;
static GLuint g_GlobalMatricesUBO;

// This gets compiled by both c and c++ compilers
#ifdef __cplusplus
extern "C" {
#endif
void init_world(world_t *state,  int screen_width, int screen_height);
void rotate_camera_yaw(world_t *statem, float degrees);
void rotate_camera_pitch(world_t *statem, float degrees);
void rotate_camera_yaw_pitch(world_t *state, float degrees_pitch, float degrees_yaw);
void zoom_view(world_t *state, float d_zoom_factor);
void move_along_view(world_t *state, float dx);
void update_view(world_t *state);
#ifdef __cplusplus
}
#endif

#endif
