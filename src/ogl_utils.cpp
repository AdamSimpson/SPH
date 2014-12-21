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

// Convert hsv to rgb
// input hsv [0:1]
// output rgb [0,1]
void hsv_to_rgb(float* HSV, float *RGB)
{
    float hue, saturation, value, hh, p, q, t, ff, r, g, b;
    long i;

    hue = HSV[0];
    saturation = HSV[1];
    value = HSV[2];

    hh = hue*360.0f;
    if(hh >= 360.0f)
        hh = 0.0f;
    hh /= 60.0f;
    i = (long)hh;
    ff = hh - i;
    p = value * (1.0f - saturation);
    q = value * (1.0f - (saturation * ff));
    t = value * (1.0f - (saturation * (1.0f - ff)));

    switch(i) {
        case 0:
            r = value;
            g = t;
            b = p;
            break;
        case 1:
            r = q;
            g = value;
            b = p;
            break;
        case 2:
            r = p;
            g = value;
            b = t;
            break;
        case 3:
            r = p;
            g = q;
            b = value;
            break;
        case 4:
            r = t;
            g = p;
            b = value;
            break;
        case 5:
            r = value;
            g = p;
            b = q;
            break;
        default:
            r = value;
            g = p;
            b = q;
            break;
    }

    RGB[0] = r;
    RGB[1] = g;
    RGB[2] = b;
}
