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

#include "glfw_utils.h"

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
void render_dividers(dividers_t *state, float *node_edges, float *colors_by_rank, int num_nodes)
{
    // 6 verticies per line, 5(x,y,r,g,b) elements per vertex, 2 lines per node
    float *verticies =  (float*)malloc(sizeof(float)*30*2*num_nodes);

    // Create verticies
    int i, offset;
    float edge_x;
    float half_width = state->screen_width*0.0000012;
    float *color;
    for(i=0; i<2*num_nodes; i++) {
        edge_x = node_edges[i];

        // Node start
        if(i%2 == 0)
            edge_x += half_width;
        else // Node end
            edge_x -= half_width;
        
        // Set color same as particle color
        color = &colors_by_rank[3*(i/2)];

        // Left triangle
        offset = i*30;
        verticies[offset]   = edge_x - half_width;
        verticies[offset+1] = -1.0f;
        verticies[offset+2] = color[0];
        verticies[offset+3] = color[1];
        verticies[offset+4] = color[2];
        verticies[offset+5] = edge_x + half_width;
        verticies[offset+6] = -1.0f;
        verticies[offset+7] = color[0];
        verticies[offset+8] = color[1];
        verticies[offset+9] = color[2];
        verticies[offset+10] = edge_x - half_width;
        verticies[offset+11] = 1.0f;
        verticies[offset+12] = color[0];
        verticies[offset+13] = color[1];
        verticies[offset+14] = color[2];

        // Right triangle
        verticies[offset+15]  = edge_x + half_width;
        verticies[offset+16]  = -1.0f;
        verticies[offset+17] = color[0];
        verticies[offset+18] = color[1];
        verticies[offset+19] = color[2];
        verticies[offset+20]  = edge_x + half_width;
        verticies[offset+21]  = 1.0f;
        verticies[offset+22] = color[0];
        verticies[offset+23] = color[1];
        verticies[offset+24] = color[2];
        verticies[offset+25] = edge_x - half_width;
        verticies[offset+26] = 1.0f;
        verticies[offset+27] = color[0];
        verticies[offset+28] = color[1];
        verticies[offset+29] = color[2];
    }

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Orphan current buffer
    glBufferData(GL_ARRAY_BUFFER, 30*2*num_nodes*sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 30*2*num_nodes*sizeof(GLfloat), verticies, GL_STREAM_DRAW);

    // Unbind buffer
    draw_dividers(state, num_nodes);

    free(verticies);
}

void create_dividers_buffers(dividers_t *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

void create_dividers_shaders(dividers_t *state)
{
    // Compile vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    compile_shader(vertexShader, "shaders/divider.vert");

    // Compile frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    compile_shader(fragmentShader, "shaders/divider.frag");

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertexShader);
    glAttachShader(state->program, fragmentShader); 

    // Link and use program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Get position location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get color location
    state->color_location = glGetAttribLocation(state->program, "color");
}

void draw_dividers(dividers_t *state, int num_nodes)
{
    // Bind circle shader program
    glUseProgram(state->program);

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->color_location, 3, GL_FLOAT, GL_FALSE, 5*sizeof(GL_FLOAT), (void*)(2*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->color_location);

    // Draw, two triangles per line
    glDrawArrays(GL_TRIANGLES, 0, 6*2*num_nodes);
}
