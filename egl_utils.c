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

/*
Copyright (c) 2012, Broadcom Europe Ltd
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the copyright holder nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "egl_utils.h"
#include "linux/input.h"
#include "renderer.h"

bool window_should_close(gl_t *state)
{
    if(state->window_should_close)
        return true;
    else
        return false;
}

// Description: Sets the display, OpenGL|ES context and screen stuff
void init_ogl(gl_t *state, render_t *render_state)
{

    bcm_host_init();

    // Initialize struct
    memset(state, 0, sizeof(gl_t));
    state->keyboard_fd = -1;
    state->mouse_fd = -1;
    state->window_should_close = false;

    // Set a user pointer up
    // This mimics GLFW
    state->user_pointer = render_state;

    int32_t success = 0;
    EGLBoolean result;
    EGLint num_config;

    static EGL_DISPMANX_WINDOW_T nativewindow;

    DISPMANX_ELEMENT_HANDLE_T dispman_element;
    DISPMANX_DISPLAY_HANDLE_T dispman_display;
    DISPMANX_UPDATE_HANDLE_T dispman_update;
    VC_RECT_T dst_rect;
    VC_RECT_T src_rect;

    static const EGLint attribute_list[] =
    {
       EGL_RED_SIZE, 8,
       EGL_GREEN_SIZE, 8,
       EGL_BLUE_SIZE, 8,
       EGL_ALPHA_SIZE, 8,
       EGL_SURFACE_TYPE, EGL_WINDOW_BIT,
       EGL_NONE
    };
   
    static const EGLint context_attributes[] =
    {
       EGL_CONTEXT_CLIENT_VERSION, 2,
       EGL_NONE
    }; 

    EGLConfig config;

    // get an EGL display connection
    state->display = eglGetDisplay(EGL_DEFAULT_DISPLAY);
    assert(state->display!=EGL_NO_DISPLAY);
 
    // initialize the EGL display connection
    result = eglInitialize(state->display, NULL, NULL);
    assert(EGL_FALSE != result);

    // get an appropriate EGL frame buffer configuration
    result = eglChooseConfig(state->display, attribute_list, &config, 1, &num_config);
//    result = eglSaneChooseConfigBRCM(state->display, attribute_list, &config, 1, &num_config);
    assert(EGL_FALSE != result);

    // get an appropriate EGL frame buffer configuration
    result = eglBindAPI(EGL_OPENGL_ES_API);
    assert(EGL_FALSE != result);

    // create an EGL rendering context
    state->context = eglCreateContext(state->display, config, EGL_NO_CONTEXT, context_attributes);
//    state->context = eglCreateContext(state->display, config, EGL_NO_CONTEXT, NULL);
    assert(state->context!=EGL_NO_CONTEXT);

    // create an EGL window surface
    success = graphics_get_display_size(0 /* LCD */, &state->screen_width, &state->screen_height);
    assert( success >= 0 );

    dst_rect.x = 0;
    dst_rect.y = 0;
    dst_rect.width = state->screen_width;
    dst_rect.height = state->screen_height;
      
    src_rect.x = 0;
    src_rect.y = 0;
    src_rect.width = state->screen_width << 16;
    src_rect.height = state->screen_height << 16;

    dispman_display = vc_dispmanx_display_open( 0 /* LCD */);
    dispman_update = vc_dispmanx_update_start( 0 );
         
    dispman_element = vc_dispmanx_element_add ( dispman_update, dispman_display,
      0/*layer*/, &dst_rect, 0/*src*/,
      &src_rect, DISPMANX_PROTECTION_NONE, 0 /*alpha*/, 0/*clamp*/, 0/*transform*/);
      
    nativewindow.element = dispman_element;
    nativewindow.width = state->screen_width;
    nativewindow.height = state->screen_height;
    vc_dispmanx_update_submit_sync( dispman_update );
      
    state->surface = eglCreateWindowSurface( state->display, config, &nativewindow, NULL );
    assert(state->surface != EGL_NO_SURFACE);

    // connect the context to the surface
    result = eglMakeCurrent(state->display, state->surface, state->surface, state->context);
    assert(EGL_FALSE != result);
}

void swap_ogl(gl_t *state)
{
    eglSwapBuffers(state->display, state->surface);
}

void exit_ogl(gl_t *state)
{
   // clear screen
   glClear( GL_COLOR_BUFFER_BIT );
   eglSwapBuffers(state->display, state->surface);

   // Release OpenGL resources
   eglMakeCurrent( state->display, EGL_NO_SURFACE, EGL_NO_SURFACE, EGL_NO_CONTEXT );
   eglDestroySurface( state->display, state->surface );
   eglDestroyContext( state->display, state->context );
   eglTerminate( state->display );

   close(state->keyboard_fd);

   printf("close\n");
}

// Return mouse position in OpenGL coordinates
// This is different than default GLFW coordinates
void get_mouse(float *x_pos, float *y_pos, gl_t *state)
{
    // Screen dimensions in pixels
    const int width = state->screen_width;
    const int height = state->screen_height;

    // Used to determine what direction dx,dy is in
    const int XSIGN = 1<<4;
    const int YSIGN = 1<<5;

    // Initialize coordinates
    static int x = 0;
    static int y = 0;

    float x_scaled;
    float y_scaled;

    ssize_t bytes_read;

    // Open file containing mouse events
    if (state->mouse_fd < 0)
        state->mouse_fd = open("/dev/input/mice", O_NONBLOCK);

    if (state->mouse_fd >= 0) {
        MOUSE_INPUT event;

        bytes_read = read(state->mouse_fd, &event, sizeof(MOUSE_INPUT));

        if(bytes_read != sizeof(MOUSE_INPUT))
            return;

        int speed_multiplier = 2;

        // Increment pixel positions
        x+=speed_multiplier*event.dx;
        y+=speed_multiplier*event.dy;

        if(event.button&XSIGN)
            x -= speed_multiplier*256;
        if(event.button&YSIGN)
            y -= speed_multiplier*256;

        // Make sure not to go out of bounds
        if(x < 0)
            x = 0.0;
        else if(x > width)
            x = width;
        if(y < 0)
            y = 0.0;
        else if(y > height)
            y = height;

        // convert to OpenGL screen coordinates from pixels
        x_scaled = (float)x/(0.5f*width) - 1.0f;
        y_scaled = (float)y/(0.5f*height) - 1.0f;

//        debug_print("pixels(%d,%d) scaled(%f,%f)\n", x,y,x_scaled,y_scaled);
    }

    *x_pos = x_scaled;
    *y_pos = y_scaled;
}


void check_key_press(gl_t *state)
{
    render_t *render_state = (render_t*)state->user_pointer;

    // Get key press
    int key = get_key_press(state);

    switch(key)
    {
        case KEY_RIGHT:
            increase_parameter(render_state);
            break;
        case KEY_LEFT:
            decrease_parameter(render_state);
            break;
        case KEY_UP:
            move_parameter_up(render_state);
            break;
        case KEY_DOWN:
            move_parameter_down(render_state);
            break;
        case KEY_PAGEUP:
            add_partition(render_state);
            break;
        case KEY_PAGEDOWN:
            remove_partition(render_state);
            break;
        case KEY_X:
            set_fluid_x(render_state);
            break;
        case KEY_Y:
            set_fluid_y(render_state);
            break;
        case KEY_A:
            set_fluid_a(render_state);
            break;
        case KEY_B:
            set_fluid_b(render_state);
            break;
        case BTN_BACK:
	    toggle_dividers(render_state);
            break;
        case BTN_FORWARD:
            toggle_dividers(render_state);
            break;
        case KEY_ESC:
            toggle_pause(render_state);
            break;
        case KEY_TAB:
            state->window_should_close = true;
            break;
 
    }
}

int get_key_press(gl_t *state)
{
    int key_code = -1;

    // Open file containing keyboard events
    // event0 or event1 depends on the USB PORT!!!!
    if (state->keyboard_fd < 0)
        state->keyboard_fd = open("/dev/input/event1",O_RDONLY|O_NONBLOCK);
    if (state->keyboard_fd >= 0) {
        struct input_event event;
        read(state->keyboard_fd, &event, sizeof(struct input_event));
	if(event.type == EV_KEY && event.value == 1 && event.code > 0)
	    key_code = (int)event.code;
    }

    return key_code;
}
