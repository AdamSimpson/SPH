#ifndef GLFW_UTILS_H
#define GLFW_UTILS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "renderer.h"

typedef struct gl_t {
    int screen_width;
    int screen_height;

    GLFWwindow* window;

} gl_t;

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void error_callback(int error, const char* description);
void check_key_press(gl_t *state);
void get_mouse(float *x, float *y, gl_t *state);
void init_ogl(gl_t *state, render_t *render_state);
void exit_ogl(gl_t *state);
void swap_ogl(gl_t *state);
bool window_should_close(gl_t *state);

#endif
