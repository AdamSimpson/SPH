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

#include "gl.hpp"
#include "renderer.hpp"
#include "tunable_parameters.hpp"

void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Get renderer object from GLFW user pointer
    Renderer *renderer = (Renderer*)glfwGetWindowUserPointer(window);

    if(action == GLFW_PRESS)
    {
        // Let renderer know of activity
        renderer->set_activity_time();

        switch(key)
        { 
            case GLFW_KEY_ESCAPE:
                self->exit_program(window);              
	        break;
            case GLFW_KEY_RIGHT:
                renderer->tunable_parameters->increase_parameter();
                break;
            case GLFW_KEY_LEFT:
                renderer->tunable_parameters->decrease_parameter();
                break;
            case GLFW_KEY_UP:
                renderer->tunable_parameters->move_parameter_up();
                break;
            case GLFW_KEY_DOWN:
                renderer->tunable_parameters->move_parameter_down();
                break;
            case GLFW_KEY_LEFT_BRACKET:
                renderer->tunable_parameters->remove_partition();
                break;
            case GLFW_KEY_RIGHT_BRACKET:
                renderer->tunable_parameters->add_partition();
                break;
            case GLFW_KEY_P:
                renderer->toggle_pause();
                break;
            case GLFW_KEY_LEFT_SHIFT:
                glfwSetCursorPos (window, renderer->tunable_parameters->cursor_view_x, renderer->tunable_parameters->cursor_view_y);
                renderer->enable_view_controls();
                break;
        }
    }
    else if (action == GLFW_RELEASE)
    {
        switch(key)
        {
            case GLFW_KEY_LEFT_SHIFT:
                glfwSetCursorPos (window, renderer->tunable_parameters->cursor_x, renderer->tunable_parameters->cursor_y);
                renderer->disable_view_controls();
                break;
        }
    }
}

static void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    // Get renderer object from GLFW user pointer
    Renderer *renderer = (Renderer*)glfwGetWindowUserPointer(window);

    // Make sure mouse stays in bounds
    // GLFW_CURSOR_DISABLED works better on mac than GLFW_CURSOR_HIDDEN but doesn't trap cursor
    // Sometimes jumping can happen when using SetCursorPos with HIDDEN mode
    float pixel_width  = renderer->tunable_parameters->screen_width;
    float pixel_height = renderer->tunable_parameters->screen_height;
    if(xpos < 0.0) {
        xpos = 0.0;
        glfwSetCursorPos(window, xpos, ypos);
    }
    else if(floor(xpos) > renderer->tunable_parameters->screen_width) {
        xpos = (double)pixel_width;
        glfwSetCursorPos(window, xpos, ypos);
    }
    if(ypos < 0.0) {
        ypos = 0.0;
        glfwSetCursorPos(window, xpos, ypos);
    }
    else if(floor(ypos) > renderer->tunable_parameters->screen_height) {
        ypos = (double)pixel_height;
        glfwSetCursorPos(window, xpos, ypos);
    }

    // Let renderer know of activity
    renderer->tunable_parameters->set_activity_time();

    if(renderer->tunable_parameters->view_controls) {
        renderer->tunable_parameters->cursor_view_x = xpos;
        renderer->tunable_parameters->cursor_view_y = ypos;

        float new_x, new_y;
        new_y = renderer->tunable_parameters->screen_height - ypos; // Flip y = 0
        new_y = new_y/(0.5*renderer->tunable_parameters->screen_height) - 1.0;
        new_x = render_state->screen_width - xpos; // Flip x = 0
        new_x = new_x/(0.5*renderer->tunable_parameters->screen_width) - 1.0;

        renderer->set_view_angle(new_x, new_y);
    }
    else if(!renderer->tunable_parameters->view_controls) {
        renderer->tunable_parameters->cursor_x = xpos;
        renderer->tunable_parameters->cursor_y = ypos;

        float new_x, new_y;
        new_y = (renderer->tunable_parameters->screen_height - ypos); // Flip y = 0
        new_y = new_y/(0.5*renderer->tunable_parameters->screen_height) - 1.0;
        new_x = xpos/(0.5*renderer->tunable_parameters->screen_width) - 1.0;

        renderer->tunable_parameters->set_mover_gl_center(new_x, new_y, -0.8f);
    }
}

// scroll wheel callback
void wheel_callback(GLFWwindow* window, double x, double y)
{
    // Get renderer object from GLFW user pointer
    Renderer *renderer = (Renderer*)glfwGetWindowUserPointer(window);

    // Let renderer know of activity
    renderer->set_activity_time();    

    if(renderer->view_controls) {
        if(x > 0.0)
            renderer->move_in_view();
        else if(x < 0.0)
            renderer->move_out_view();
    }
    else {
        // Call increase/decrease mover calls
        if(y > 0.0)
	    renderer->tunable_parameters->increase_mover_radius();
        else if(y < 0.0)
	    renderer->tunable_parameters->decrease_mover_radius();
        if(x > 0.0)
            renderer->tunable_parameters->increase_mover_radius();
        else if(x < 0.0)
            renderer->tunable_parameters->decrease_mover_radius();
    }
}

// Description: Sets the display, OpenGL context and screen stuff
void GL::GL()
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

    glfwGetFramebufferSize(this->window, this->screen_width, this->screen_height);
    glViewport(0, 0, this->screen_width, this->screen_height);

    if(!state->window)
	exit(EXIT_FAILURE);

    // Set current context to window
    glfwMakeContextCurrent(this->window);

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
    glfwSetKeyCallback(this->window, key_callback);

    // Set mouse cursor callback
    glfwSetCursorPosCallback(this->window, mouse_callback);

    // Set scroll wheel callback
    glfwSetScrollCallback(this->window, wheel_callback);

    // Add render state to window pointer
    // Used for key callbacks
    glfwSetWindowUserPointer(this->window, render_state);

    // Disable regular cursor
    glfwSetInputMode(this->window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // Disable vsync for true FPS testing
    // Default limit 60 fps
//    glfwSwapInterval(0);

    // Set background color and clear buffers
    glClearColor(0.15f, 0.25f, 0.35f, 1.0f);
    glClear( GL_COLOR_BUFFER_BIT );
}

bool GL::window_should_close()
{
    if(glfwWindowShouldClose(self->window))
        return true;
    else
        return false;
}

void GL::check_user_input()
{
    // Poll GLFW for key press or mouse input
    glfwPollEvents();
}

void GL::swap_buffers()
{
    glfwSwapBuffers(this->window);
}

void GL::exit()
{
//    glDeleteProgram(state->shaderProgram);
//    glDeleteShader(fragmentShader);
//    glDeleteShader(vertexShader);
//    glDeleteBuffers(1, &vbo);
//    glDeleteVertexArrays(1, &vao);

    glfwDestroyWindow(this->window);
    glfwTerminate();

    printf("close\n");
}

// Convert pixel coordinates, lower left origin, to gl coordinates, center origin
void GL::pixel_to_gl(int pixel_x, int pixel_y, float *gl_x, float *gl_y)
{
    float half_x = this->screen_width/2.0;
    float half_y = this->screen_height/2.0;
    *gl_x = pixel_x/half_x - 1.0;
    *gl_y = pixel_y/half_y - 1.0;

}

// Exit and set return value for specific program if one selected
void GL::exit_program(GLFWwindow* window)
{
    glfwSetWindowShouldClose(window, GL_TRUE);
}

