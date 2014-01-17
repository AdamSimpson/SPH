#include <stdio.h>
#include <assert.h>

#include "egl_utils.h"

#include "GLES2/gl2.h"
#include "EGL/egl.h"
#include "EGL/eglext.h"

#include "circles_gl.h"

#include "bcm_host.h"

void update_points(float *points, int num_points, STATE_T *state)
{
    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 5*num_points*sizeof(GLfloat), points, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    draw_circles(state, num_points);
}

void create_shaders(STATE_T *state)
{
    // Shader source
    const GLchar* vertexSource =
        "attribute vec2 position;"
        "attribute vec3 color;"
	"varying vec4 frag_color;"
        "void main() {"
        "   gl_Position = vec4(position, 0.0, 1.0);"
        "   gl_PointSize = 10.0;"
	"   frag_color = vec4(color, 1.0);"
        "}";
    const GLchar* fragmentSource =
        "precision mediump float;"
	"varying vec4 frag_color;"
        "const mediump vec2 center = vec2(0.5, 0.5);"
        "const mediump float radius = 0.5;"
        "void main() {"
	"    mediump float distance_from_center = distance(center, gl_PointCoord);"
	"    lowp float in_radius = step(distance_from_center, radius);"
        "    gl_FragColor = frag_color * in_radius;"
        "}";

    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, NULL);
    glCompileShader(vertexShader);

    showlog(vertexShader);   

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, NULL);
    glCompileShader(fragmentShader);

    showlog(fragmentShader);

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertexShader);
    glAttachShader(state->program, fragmentShader); 
  
    // Link and use program
    glLinkProgram(state->program);
    glUseProgram(state->program);
    check();

    // Get position location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get tex_coord location
    state->color_location = glGetAttribLocation(state->program, "color");

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


}

void draw_circles(STATE_T *state, int num_points)
{
    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->color_location, 3, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT),(void*)(2*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->color_location);

    // Draw
    glDrawArrays(GL_POINTS, 0, num_points);
}
