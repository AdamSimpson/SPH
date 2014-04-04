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
#include <stdlib.h>
#include "fluid.h"
#include "renderer.h"
#include "lodepng.h"
#include "background_gl.h"

#include "ogl_utils.h"

void create_backround_program(background_t *state)
{
    // Compile vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
        compile_shader(vertex_shader, "SPH/shaders/background_es.vert");
    #else
        compile_shader(vertex_shader, "shaders/background.vert");
    #endif

    // Compile fragment shader
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
        compile_shader(frag_shader, "SPH/shaders/background_es.frag");
    #else
        compile_shader(frag_shader, "shaders/background.frag");
    #endif

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertex_shader);
    glAttachShader(state->program, frag_shader);

    // Link  program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Get position attribute location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get tex_coord location
    state->tex_coord_location = glGetAttribLocation(state->program, "tex_coord");
    // Get tex uniform location
    state->tex_location = glGetUniformLocation(state->program, "tex");
}

void create_background_buffers(background_t *state)
{
    // VAO is required for OpenGL 3+ when using VBO I believe
    #ifndef RASPI
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);

    // Generate element buffer
    glGenBuffers(1, &state->ebo);
}

void create_background_vertices(background_t *state)
{
    // Vertices: Pos(x,y) Tex(x,y)
    // For simplicity only single vbo is generated and offset used as needed

    // Background image dimensions in gl screen coordinates
    float bg_width = 2.0*(state->background_width/(float)state->screen_width);
    float bg_height =  2.0*(state->background_height/(float)state->screen_height);

    float lower_left_x = -bg_width/2.0;
    float lower_left_y = -bg_height/2.0;
    float lower_right_x = lower_left_x + bg_width;
    float lower_right_y = lower_left_y;
    float upper_right_x = lower_right_x;
    float upper_right_y = lower_right_y + bg_height;
    float upper_left_x = lower_left_x;
    float upper_left_y = lower_left_y + bg_height;

    float vertices[] = {
         // Full screen vertices
        upper_left_x,  upper_left_y, 0.0f, 0.0f, // Upper left
        upper_right_x,  upper_right_y, 1.0f, 0.0f, // Upper right
        lower_right_x, lower_right_y, 1.0f, 1.0f, // Lower right
	lower_left_x, lower_left_y, 0.0f, 1.0f  // Lower left
    };

    printf("Screen resolution %d, %d\n", state->screen_width, state->screen_height);
    printf("ul (%f,%f), ur (%f, %f), ll (%f, %f), lr (%f, %f)\n", upper_left_x, upper_left_y, upper_right_x, upper_right_y, lower_left_x, lower_left_y, lower_right_x, lower_right_y);

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 3*4*4*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

    // Elements
    GLubyte elements[] = {
        2, 3, 0,
        0, 1, 2
    };

    // Set buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state->ebo);
    // Fill buffer
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2*3*sizeof(GLubyte), elements, GL_STATIC_DRAW);
}

void create_background_texture(background_t *state)
{
    // Read in background PNG
    unsigned error;
    unsigned char* image;
    unsigned width, height;

    #ifdef RASPI
    error = lodepng_decode32_file(&image, &width, &height, "SPH/OakRidgeLeaf.png");
    #else
    error = lodepng_decode32_file(&image, &width, &height, "OakRidgeLeaf.png");
    #endif
    if(error) printf("error %u: %s\n", error, lodepng_error_text(error));

    state->background_width = width;
    state->background_height = height;

    printf("Background image loaded: %d x %d pixels\n", width, height);

    // Generate texture
    glGenTextures(1, &state->tex_uniform);

    // Set texture unit 0 and bind texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);

    // Buffer texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);

    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Create Mipmap
    glGenerateMipmap(GL_TEXTURE_2D);

    // Release image host memory
    free(image);
}

void init_background(background_t *state, int screen_width, int screen_height)
{
    // Set screen with/height in pixels
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Create program
    create_backround_program(state);

    // Generate buffers 
    create_background_buffers(state);  

    // Create texture from image
    // Must be called before create_background_verticies
    // so image dimensions known
    create_background_texture(state);

    // Set verticies
    create_background_vertices(state);
}

void draw_background(background_t *state)
{
    // Setup program
    glUseProgram(state->program);

    // Setup buffers
    size_t vert_size = 4*sizeof(GL_FLOAT);
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, vert_size, 0);
    glEnableVertexAttribArray(state->position_location);
    glVertexAttribPointer(state->tex_coord_location, 2, GL_FLOAT, GL_FALSE, vert_size,(void*)(2*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->tex_coord_location);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state->ebo);

    // Disable Blend
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Setup texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->tex_location, 0);

    // Draw background
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);
}
