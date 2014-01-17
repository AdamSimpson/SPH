#ifndef EGL_UTILS_H
#define EGL_UTILS_H

#include "GLES2/gl2.h"
#include "EGL/egl.h"
#include "EGL/eglext.h"

typedef struct {
    uint32_t screen_width;
    uint32_t screen_height;

    EGLDisplay display;
    EGLSurface surface;
    EGLContext context;

    int keyboard_fd;
} EGL_STATE_T;


void init_ogl(EGL_STATE_T *state);
void exit_ogl(EGL_STATE_T *state);
void showlog(GLint shader);
void egl_swap(EGL_STATE_T *state);
void check();
int get_key_press(EGL_STATE_T *state);

#endif
