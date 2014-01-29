#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "glfw_utils.h"

void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

// Return mouse position in OpenGL screen coordinates
// x,y [-1, 1], center of screen is origin
void get_mouse(double *x, double *y, GL_STATE_T *state)
{
    glfwGetCursorPos(state->window, x, y);
    *y = (state->screen_height - *y); // Flip y = 0
    *y = *y/(0.5*state->screen_height) - 1.0;
    *x = *x/(0.5*state->screen_width) - 1.0;
}

// Description: Sets the display, OpenGL context and screen stuff
void init_ogl(GL_STATE_T *state)
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

    state->screen_width = 800;
    state->screen_height = 800;
    state->window = glfwCreateWindow(state->screen_width, state->screen_height, "SPH", NULL, NULL);
    if(!state->window)
	exit(EXIT_FAILURE);

    // Set current context to window
    glfwMakeContextCurrent(state->window);
    // Set key callback
    glfwSetKeyCallback(state->window, key_callback);

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    glewInit();

    // Set background color and clear buffers
    glClearColor(0.15f, 0.25f, 0.35f, 1.0f);
    glClear( GL_COLOR_BUFFER_BIT );
}

void swap_ogl(GL_STATE_T *state)
{
    glfwSwapBuffers(state->window);

    glfwPollEvents();
}

void exit_ogl(GL_STATE_T *state)
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
