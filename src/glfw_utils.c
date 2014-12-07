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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "glfw_utils.h"
#include "renderer.h"
#include "controls.h"

void check_user_input(gl_t *state)
{
    // Poll GLFW for key press or mouse input
    glfwPollEvents();
}

void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}


bool window_should_close(gl_t *state)
{
    if(glfwWindowShouldClose(state->window))
	return true;
    else
	return false;
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Get render_state from GLFW user pointer
    render_t *render_state = (render_t*)glfwGetWindowUserPointer(window);

    if(action == GLFW_PRESS)
    {
        // Let renderer know of activity
        set_activity_time(render_state);

        switch(key)
        { 
            case GLFW_KEY_ESCAPE:
                exit_program(window);              
	        break;
            case GLFW_KEY_RIGHT:
                increase_parameter(render_state);
                break;
            case GLFW_KEY_LEFT:
                decrease_parameter(render_state);
                break;
            case GLFW_KEY_UP:
                move_parameter_up(render_state);
                break;
            case GLFW_KEY_DOWN:
                move_parameter_down(render_state);
                break;
            case GLFW_KEY_LEFT_BRACKET:
                remove_partition(render_state);
                break;
            case GLFW_KEY_RIGHT_BRACKET:
                add_partition(render_state);
                break;
            case GLFW_KEY_P:
                toggle_pause(render_state);
                break;
            case GLFW_KEY_LEFT_SHIFT:
                glfwSetCursorPos (window, render_state->gl_state->cursor_view_x, render_state->gl_state->cursor_view_y);
                enable_view_controls(render_state);
                break;
        }
    }
    else if (action == GLFW_RELEASE)
    {
        switch(key)
        {
            case GLFW_KEY_LEFT_SHIFT:
                glfwSetCursorPos (window, render_state->gl_state->cursor_x, render_state->gl_state->cursor_y);
                disable_view_controls(render_state);
                break;
        }
    }
}

static void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    printf("xpos: %f, ypos: %f\n", xpos, ypos);
    // Get render_state from GLFW user pointer
    render_t *render_state = (render_t*)glfwGetWindowUserPointer(window);

    // Make sure mouse stays in bounds
    // GLFW_CURSOR_DISABLED works better on mac than GLFW_CURSOR_HIDDEN but doesn't trap cursor
    // Sometimes jumping can happen when using SetCursorPos with HIDDEN mode
    float pixel_width = render_state->gl_state->screen_width;
    float pixel_height = render_state->gl_state->screen_height;
    if(xpos < 0.0) {
        xpos = 0.0;
        glfwSetCursorPos(window, xpos, ypos);
    }
    else if(floor(xpos) > render_state->gl_state->screen_width) {
        xpos = (double)pixel_width;
        glfwSetCursorPos(window, xpos, ypos);
    }
    if(ypos < 0.0) {
        ypos = 0.0;
        glfwSetCursorPos(window, xpos, ypos);
    }
    else if(floor(ypos) > render_state->gl_state->screen_height) {
        ypos = (double)pixel_height;
        glfwSetCursorPos(window, xpos, ypos);
    }

    // Let renderer know of activity
    set_activity_time(render_state);

    if(render_state->view_controls) {
        render_state->gl_state->cursor_view_x = xpos;
        render_state->gl_state->cursor_view_y = ypos;

        float new_x, new_y;
        new_y = render_state->screen_height - ypos; // Flip y = 0
        new_y = new_y/(0.5*render_state->screen_height) - 1.0;
        new_x = render_state->screen_width - xpos; // Flip x = 0
        new_x = new_x/(0.5*render_state->screen_width) - 1.0;

        set_view_angle(render_state, new_x, new_y);
    }
    else if(!render_state->view_controls) {
        printf("move: %f, %f\n", xpos, ypos);
        render_state->gl_state->cursor_x = xpos;
        render_state->gl_state->cursor_y = ypos;

        float new_x, new_y;
        new_y = (render_state->screen_height - ypos); // Flip y = 0
        new_y = new_y/(0.5*render_state->screen_height) - 1.0;
        new_x = xpos/(0.5*render_state->screen_width) - 1.0;

        set_mover_gl_center(render_state, new_x, new_y, -0.8f);
    }
}

// scroll wheel callback
void wheel_callback(GLFWwindow* window, double x, double y)
{
    // Get render_state from GLFW user pointer
    render_t *render_state = (render_t*)glfwGetWindowUserPointer(window);

    // Let renderer know of activity
    set_activity_time(render_state);    

    if(render_state->view_controls) {
        if(x > 0.0)
            zoom_in_view(render_state);
        else if(x < 0.0)
            zoom_out_view(render_state);
    }
    else {
        // Call increase/decrease mover calls
        if(y > 0.0)
	    increase_mover_height(render_state);
        else if(y < 0.0)
	    decrease_mover_height(render_state);
        if(x > 0.0)
            increase_mover_width(render_state);
        else if(x < 0.0)
            decrease_mover_width(render_state);
    }
}

// Description: Sets the display, OpenGL context and screen stuff
void init_ogl(gl_t *state, render_t *render_state)
{
    // Set error callback
    glfwSetErrorCallback(error_callback);

    // Initialize GLFW
    if(!glfwInit())
        exit(EXIT_FAILURE);

    // Create window
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    // Hit at multisample
    glfwWindowHint(GLFW_SAMPLES, 4);
    printf("Enabling multisample! \n");

    // Enable SRGB for gamma correction
    printf("Enabling SRGB framebuffer! \n");
    glfwWindowHint(GLFW_SRGB_CAPABLE, GL_TRUE);

    // Try to get full screen
    // Retina screen is a pain...
    GLFWvidmode *mode = (GLFWvidmode*)glfwGetVideoMode(glfwGetPrimaryMonitor());
    state->window = glfwCreateWindow(mode->width, mode->height, "SPH", glfwGetPrimaryMonitor(), NULL);

    glfwGetFramebufferSize(state->window, &state->screen_width, &state->screen_height);
    glViewport(0, 0, state->screen_width, state->screen_height);

    if(!state->window)
	exit(EXIT_FAILURE);

    // Set current context to window
    glfwMakeContextCurrent(state->window);

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    glewInit();

    // Enable depth testing
    glEnable(GL_DEPTH_TEST);

    // Enable multisampling
    glEnable( GL_MULTISAMPLE );

    // Enable SRGB framebuffer
    glEnable(GL_FRAMEBUFFER_SRGB);

    // Set key callback
    glfwSetKeyCallback(state->window, key_callback);

    // Set mouse cursor callback
    glfwSetCursorPosCallback(state->window, mouse_callback);

    // Set scroll wheel callback
    glfwSetScrollCallback(state->window, wheel_callback);

    // Add render state to window pointer
    // Used for key callbacks
    glfwSetWindowUserPointer(state->window, render_state);

    // Disable regular cursor
    glfwSetInputMode(state->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

//    glfwSetInputMode(state->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Disable vsync for true FPS testing
    // Default limit 60 fps
//    glfwSwapInterval(0);

    // Set background color and clear buffers
    glClearColor(0.15f, 0.25f, 0.35f, 1.0f);
    glClear( GL_COLOR_BUFFER_BIT );
}

void swap_ogl(gl_t *state)
{
    glfwSwapBuffers(state->window);

//    glfwPollEvents();
}

void exit_ogl(gl_t *state)
{
//    glDeleteProgram(state->shaderProgram);
//    glDeleteShader(fragmentShader);
//    glDeleteShader(vertexShader);
//    glDeleteBuffers(1, &vbo);
//    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(state->window);
    glfwTerminate();

    printf("close\n");
}

// Convert pixel coordinates, lower left origin, to gl coordinates, center origin
void pixel_to_gl(gl_t *state, int pixel_x, int pixel_y, float *gl_x, float *gl_y)
{
    float half_x = state->screen_width/2.0;
    float half_y = state->screen_height/2.0;
    *gl_x = pixel_x/half_x - 1.0;
    *gl_y = pixel_y/half_y - 1.0;

}

// Exit and set return value for specific program if one selected
void exit_program(GLFWwindow* window)
{
    glfwSetWindowShouldClose(window, GL_TRUE);
}

