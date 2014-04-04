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
#include "dividers_gl.h"
#include "ogl_utils.h"

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

void init_dividers(dividers_t *state, int screen_width, int screen_height)
{
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Create circle buffers
    create_dividers_buffers(state);

    // Create and set dividers shaders
    // Also links dividers program
    create_dividers_shaders(state);
}

// Update coordinate of dividers and render
void render_dividers(dividers_t *state, float *node_edges, int num_nodes)
{
    // 6 verticies per line, 5(x,y,r,g,b) elements per vertex, 2 lines per node
    float *verticies =  (float*)malloc(sizeof(float)*2*6*2*num_nodes);

    // Create verticies
    int i, offset;
    float edge_x;
    for(i=0; i<2*num_nodes; i++) {
        float edge = node_edges[i];
        float half_width = 0.001;

        // Left triangle
        offset = i*12;
        verticies[offset]   = edge_x - half_width;
        verticies[offset+1] = -1.0f;
        verticies[offset+2] = edge_x + half_width;
        verticies[offset+3] = -1.0f;
        verticies[offset+4] = edge_x - half_width;
        verticies[offset+5] = 1.0f;

        // Right triangle
        verticies[offset+6]  = edge_x + half_width;
        verticies[offset+6]  = -1.0f;
        verticies[offset+7]  = edge_x + half_width;
        verticies[offset+8]  = 1.0f;
        verticies[offset+10] = edge_x - half_width;
        verticies[offset+11] = 1.0f;
    }

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Orphan current buffer
    glBufferData(GL_ARRAY_BUFFER, 6*2*2*num_nodes*sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 6*2*2*num_nodes*sizeof(GLfloat), verticies, GL_STREAM_DRAW);

    // Unbind buffer
    draw_dividers(state, num_nodes);
}

void create_dividers_buffers(dividers_t *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    #ifndef RASPI
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

void create_dividers_shaders(dividers_t *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
      compile_shader(vertexShader, "SPH/shaders/dividers_es.vert");
    #else
      compile_shader(vertexShader, "shaders/dividers.vert");
    #endif

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
      compile_shader(fragmentShader, "SPH/shaders/dividers_es.frag");
    #else
      compile_shader(fragmentShader, "shaders/dividers.frag");
    #endif

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertexShader);
    glAttachShader(state->program, fragmentShader); 

    // Link and use program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Get position location
    state->position_location = glGetAttribLocation(state->program, "position");
}

void draw_dividers(dividers_t *state, int num_nodes)
{
    // Bind circle shader program
    glUseProgram(state->program);

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw, two triangles per line
    glDrawArrays(GL_TRIANGLES, 0, 2*2*num_nodes);
}
