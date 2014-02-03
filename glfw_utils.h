#ifndef GLFW_UTILS_H
#define GLFW_UTILS_H

#include <GL/glew.h>
#include <GLFW/glfw3.h>


typedef struct {
    uint32_t screen_width;
    uint32_t screen_height;

    GLFWwindow* window;

} GL_STATE_T;

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void error_callback(int error, const char* description);
void check_key_press(GL_STATE_T *state);
void get_mouse(double *x, double *y, GL_STATE_T *state);
void init_ogl(GL_STATE_T *state);
void exit_ogl(GL_STATE_T *state);
void swap_ogl(GL_STATE_T *state);

#endif
