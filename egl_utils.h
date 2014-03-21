#ifndef EGL_UTILS_H
#define EGL_UTILS_H

#include "GLES2/gl2.h"
#include "EGL/egl.h"
#include "EGL/eglext.h"

#include "fcntl.h"
#include "linux/input.h"

#include "bcm_host.h"
#include "renderer.h"


typedef struct gl_t {
    uint32_t screen_width;
    uint32_t screen_height;

    EGLDisplay display;
    EGLSurface surface;
    EGLContext context;

    int keyboard_fd;
    int mouse_fd;

    void *user_pointer; // mimics GLFW user pointer

    bool window_should_close;
} gl_t;

typedef struct {
    char button;
    char dx;
    char dy;
} MOUSE_INPUT;

void init_ogl(gl_t *state, RENDER_T *render_state);
void exit_ogl(gl_t *state);
void swap_ogl(gl_t *state);
void check_key_press(gl_t *state);
int get_key_press(gl_t *state);
void get_mouse(float *x_pos, float *y_pos, gl_t *state);
bool window_should_close(gl_t *state);

#endif
