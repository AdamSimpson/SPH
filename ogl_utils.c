#include <stdio.h>
#include <assert.h>
#include <sys/stat.h>
#include <stdlib.h>
#include "ogl_utils.h"

inline void check()
{
    GLenum err = glGetError();
    if(err != GL_NO_ERROR) {
        printf("GL Error: %d\n", err);
        exit(EXIT_FAILURE);
    }
}

void show_program_log(GLuint program)
{
    char log[1024];
    glGetProgramInfoLog(program, sizeof(log), NULL, log);
    printf("%d:program:\n%s\n", program, log);

}

void showlog(GLint shader)
{
   // Prints the compile log for a shader
   char log[1024];
   glGetShaderInfoLog(shader,sizeof(log), NULL, log);
   printf("%d:shader:\n%s\n", shader, log);
}

void compile_shader(GLuint shader, const char *file_name)
{
    // Read shader source from file_name
    FILE *fh = fopen(file_name, "r");
    if(!fh) {
        printf("Error: Failed to open shader\n");
    }
    struct stat statbuf;
    stat(file_name, &statbuf);
    char *shader_source = (char *) malloc(statbuf.st_size + 1);
    fread(shader_source, statbuf.st_size, 1, fh);
    shader_source[statbuf.st_size] = '\0';

    // Compile shader
    const GLchar *gl_shader_source = shader_source;
    glShaderSource(shader, 1, &gl_shader_source, NULL);
    glCompileShader(shader);

    showlog(shader);

    free(shader_source);
}

