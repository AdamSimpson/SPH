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

#ifndef EGL_UTILS_H
#define EGL_UTILS_H

#include "GLES2/gl2.h"
#include "EGL/egl.h"
#include "EGL/eglext.h"

#include "fcntl.h"
#include "linux/input.h"

#include "bcm_host.h"
#include "renderer.h"
#include "controls.h"
#include <stdbool.h>

typedef struct gl_t {
    uint32_t screen_width;
    uint32_t screen_height;

    EGLDisplay display;
    EGLSurface surface;
    EGLContext context;

    int controller_1_fd;
    int controller_2_fd;

    void *user_pointer; // mimics GLFW user pointer

    bool window_should_close;
} gl_t;

typedef struct {
    char button;
    char dx;
    char dy;
} MOUSE_INPUT;

#ifdef LEAP_MOTION_ENABLED
  #include <curl/curl.h>

  typedef struct curl_mem_t {
    char *memory;
    size_t size;
  } curl_mem_t;

  CURL *curl_handle;
  CURLcode curl_res;
  curl_mem_t curl-chunk;

  static size_t curl_callback(void *contents, size_t size, size_t count, void *user_p);
#endif


void init_ogl(gl_t *state, render_t *render_state);
void process_controller_events(gl_t *state, int controller_fd);
void exit_ogl(gl_t *state);
void swap_ogl(gl_t *state);
void check_user_input(gl_t *state);
void handle_key(gl_t *state, struct input_event *event);
void handle_mouse(gl_t *state, struct input_event *event);
void handle_joystick(gl_t *state, struct input_event *event);
bool window_should_close(gl_t *state);
void pixel_to_gl(gl_t *state, int pixel_x, int pixel_y, float *gl_x, float *gl_y);
void exit_with_selected_program(render_t *render_state, gl_t *gl_state);

#endif
