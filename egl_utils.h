#ifndef EGL_UTILS_H
#define EGL_UTILS_H

#include "GLES2/gl2.h"
#include "EGL/egl.h"
#include "EGL/eglext.h"

#include "fcntl.h"
#include "linux/input.h"

#include "bcm_host.h"

typedef struct {
    uint32_t screen_width;
    uint32_t screen_height;

    EGLDisplay display;
    EGLSurface surface;
    EGLContext context;

    int keyboard_fd;
    int mouse_fd;
} GL_STATE_T;

typedef struct {
    char button;
    char dx;
    char dy;
} MOUSE_INPUT;

void init_ogl(GL_STATE_T *state);
void exit_ogl(GL_STATE_T *state);
void swap_ogl(GL_STATE_T *state);
void check_key_press(GL_STATE_T *state);
int get_key_press(GL_STATE_T *state);
void get_mouse(double *x_pos, double *y_pos, GL_STATE_T *state);

#endif
