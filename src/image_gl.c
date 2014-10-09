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
#include "lodepng.h"
#include "image_gl.h"

#include "ogl_utils.h"

void create_image_program(image_t *state)
{
    // Compile vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
        compile_shader(vertex_shader, "SPH/shaders/image_es.vert");
    #else
        compile_shader(vertex_shader, "shaders/image.vert");
    #endif

    // Compile fragment shader
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
        compile_shader(frag_shader, "SPH/shaders/image_es.frag");
    #else
        compile_shader(frag_shader, "shaders/image.frag");
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

void create_image_buffers(image_t *state)
{
    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);

    // Generate element buffer
    glGenBuffers(1, &state->ebo);
}

void create_image_vertices(image_t *state)
{
    // Vertices: Pos(x,y) Tex(x,y)
    // For simplicity only single vbo is generated and offset used as needed

    // image dimensions in gl screen coordinates
    float image_width = 2.0*(state->image_width/(float)state->gl_state->screen_width);
    float image_height =  2.0*(state->image_height/(float)state->gl_state->screen_height);

    float gl_x, gl_y;
//    pixel_to_gl(state->gl_state, state->lower_left_x, state->lower_left_y, &gl_x, &gl_y);

    gl_x = state->lower_left_x;
    gl_y = state->lower_left_y;

    float lower_left_x = gl_x;
    float lower_left_y = gl_y;
    float lower_right_x = lower_left_x + image_width;
    float lower_right_y = lower_left_y;
    float upper_right_x = lower_right_x;
    float upper_right_y = lower_right_y + image_height;
    float upper_left_x = lower_left_x;
    float upper_left_y = lower_left_y + image_height;

    float vertices[] = {
         // Full screen vertices
        upper_left_x,  upper_left_y, 0.0f, 0.0f, // Upper left
        upper_right_x,  upper_right_y, 1.0f, 0.0f, // Upper right
        lower_right_x, lower_right_y, 1.0f, 1.0f, // Lower right
	lower_left_x, lower_left_y, 0.0f, 1.0f  // Lower left
    };

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

void create_image_texture(image_t *state)
{
    // Read in image PNG
    unsigned error;
    unsigned char* image;
    unsigned width, height;

    error = lodepng_decode32_file(&image, &width, &height, state->file_name);
    if(error) printf("error %u: %s\n", error, lodepng_error_text(error));

    printf("loaded image: %dx%d\n", width, height);

    // Photoshop sets transparent pixels to (1,1,1,0) in PNG
    // When averaged mipmaps are generated this will produce white border artifacts
    // Setting any 0 transparency pixels to RGB = 0 will fix this bleeding
    int i;
    for(i=3; i<height*width*4; i+=4){
        if(image[i] == 0) {
            image[i-1] = 0;
            image[i-2] = 0;
            image[i-3] = 0;
        }
    }

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

    // Unbind texture
    glBindTexture(GL_TEXTURE_2D, 0);

    error = lodepng_decode32_file(&image, &width, &height, state->selected_file_name);
    if(error) printf("error %u: %s\n", error, lodepng_error_text(error));

    printf("loaded image: %dx%d\n", width, height);

    // Photoshop sets transparent pixels to (1,1,1,0) in PNG
    // When averaged mipmaps are generated this will produce white border artifacts
    // Setting any 0 transparency pixels to RGB = 0 will fix this bleeding
    for(i=3; i<height*width*4; i+=4){
        if(image[i] == 0) {
            image[i-1] = 0;
            image[i-2] = 0;
            image[i-3] = 0;
        }
    }

    // Generate texture
    glGenTextures(1, &state->tex_selected_uniform);

    // Set texture unit 0 and bind texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_selected_uniform);

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

    // Unbind texture
    glBindTexture(GL_TEXTURE_2D, 0);
}


void init_image(image_t *state, gl_t *gl_state, char *file_name, char *selected_file_name, float lower_left_x, float lower_left_y, int image_width, int image_height)
{
    // Set GL state
    state->gl_state = gl_state;

    // Set selected to false
    state->selected = false;

    // Set filename
    strcpy(state->file_name, file_name);
    strcpy(state->selected_file_name, selected_file_name);

    // Set lower left image pixel position
    state->lower_left_x = lower_left_x;
    state->lower_left_y = lower_left_y;

    // Image display dimensions in pixels
    state->image_width = image_width;
    state->image_height = image_height;

    // Create program
    create_image_program(state);

    // Generate buffers 
    create_image_buffers(state);  

    // Create texture from image
    // Must be called before create_image_verticies
    // so image dimensions known
    create_image_texture(state);

    // Set verticies
    create_image_vertices(state);
}

void draw_image(image_t *state)
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

    // Enable Blend
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Setup texture
    glActiveTexture(GL_TEXTURE0);
    if(state->selected)
        glBindTexture(GL_TEXTURE_2D, state->tex_selected_uniform);
    else
        glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->tex_location, 0);

    // Draw image
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

    // Unbind texture
    glBindTexture(GL_TEXTURE_2D, 0);

    glDisable(GL_BLEND);
}

// test if x,y gl position is inside of image
bool inside_image(image_t *state, float gl_x, float gl_y)
{
    // image dimensions in gl screen coordinates
    float image_width = 2.0f*(state->image_width/(float)state->gl_state->screen_width);
    float image_height =  2.0f*(state->image_height/(float)state->gl_state->screen_height);

    if(gl_x > state->lower_left_x && gl_x < (state->lower_left_x + image_width) ){
        if(gl_y > state->lower_left_y && gl_y < (state->lower_left_y + image_height) ) {
            return true;
        }
    }

    return false;
}
