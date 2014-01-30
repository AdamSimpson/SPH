#ifndef OGL_UTILS_H
#define OGL_UTILS_H

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

inline void check();
void showlog(GLint shader);
void compile_shader(GLuint shader, const char *file_name);

#endif
