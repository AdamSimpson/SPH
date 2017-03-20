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
#include "exit_menu_gl.h"

static float screen_scale = 1;

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
    render_t *render_state = glfwGetWindowUserPointer(window);

    if(action == GLFW_PRESS || action == GLFW_REPEAT)
    {

        // Let renderer know of activity
        set_activity_time(render_state);

        switch(key)
        { 
            case GLFW_KEY_ESCAPE:
                #ifdef EXIT_SIMPLE
                exit_now(render_state, window);
                #else
                toggle_quit_mode(render_state);
                #endif        
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
            case GLFW_KEY_X:
                set_fluid_x(render_state);
                break;
            case GLFW_KEY_Y:
                set_fluid_y(render_state);
                break;
            case GLFW_KEY_A:
                if(render_state->quit_mode)
                    exit_with_selected_program(render_state, window);
                set_fluid_a(render_state);
                break;
            case GLFW_KEY_B:
                set_fluid_b(render_state);
                break;
            case GLFW_KEY_D:
                toggle_dividers(render_state);
                break;
            case GLFW_KEY_P:
                toggle_pause(render_state);
                break;
            case GLFW_KEY_L:
                toggle_liquid(render_state);
                break;
        }
    }
}

static void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    // Get render_state from GLFW user pointer
    render_t *render_state = glfwGetWindowUserPointer(window);

    // Let renderer know of activity
    set_activity_time(render_state);

    float new_x, new_y;
    new_y = (render_state->screen_height*screen_scale - ypos); // Flip y = 0
    new_y = new_y/(0.5*render_state->screen_height*screen_scale) - 1.0;
    new_x = xpos/(0.5*render_state->screen_width*screen_scale) - 1.0;

    set_mover_gl_center(render_state, new_x, new_y);
}

// scroll wheel callback
void wheel_callback(GLFWwindow* window, double x, double y)
{
    // Get render_state from GLFW user pointer
    render_t *render_state = glfwGetWindowUserPointer(window);

    // Let renderer know of activity
    set_activity_time(render_state);    

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

    // Try to get full screen
    // Retina screen is a pain...
    GLFWvidmode *mode = (GLFWvidmode*)glfwGetVideoMode(glfwGetPrimaryMonitor());
    state->window = glfwCreateWindow(mode->width, mode->height, "SPH", glfwGetPrimaryMonitor(), NULL);

    int window_width, window_height;
    glfwGetWindowSize(state->window, &window_width, &window_height);
    glfwGetFramebufferSize(state->window, &state->screen_width, &state->screen_height);
    glViewport(0, 0, state->screen_width, state->screen_height);

    // Retina screens don't have a 1 to 1 ratio 
    screen_scale = window_width / (float)state->screen_width;

    if(!state->window)
	exit(EXIT_FAILURE);

    // Set current context to window
    glfwMakeContextCurrent(state->window);

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    glewInit();

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
    glfwSetInputMode(state->window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

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
void exit_with_selected_program(render_t *render_state, GLFWwindow* window)
{
    if(render_state->exit_menu_state->mandelbrot_state->selected) {
        glfwSetWindowShouldClose(window, GL_TRUE);
        render_state->return_value = 10;
    }
    else if (render_state->exit_menu_state->sph_state->selected) {
        glfwSetWindowShouldClose(window, GL_TRUE);
        render_state->return_value = 20;
    }
    else if (render_state->exit_menu_state->terminal_state->selected) {
        glfwSetWindowShouldClose(window, GL_TRUE);
        render_state->return_value = 0;
    }
}

// Exit without using program selector
void exit_now(render_t *render_state, GLFWwindow* window) {
  glfwSetWindowShouldClose(window, GL_TRUE);
  render_state->return_value = 0;
}
