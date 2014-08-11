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
#include "liquid_gl.h"
#include "ogl_utils.h"

#ifdef GLFW
  #include "glfw_utils.h"
#else
  #include "egl_utils.h"
#endif

void init_liquid(liquid_t *state, int screen_width, int screen_height)
{
    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Amount fluid texture will be reduced from screen resolution
    #ifdef RASPI
    state->reduction = 8;
    #else
    state->reduction = 2;
    #endif

    // Create circle buffers
    create_liquid_buffers(state);

    // Create and set particle shaders
    // Also links particle program
    create_liquid_shaders(state);

    // Set verticies
    create_texture_verticies(state);
}

// Update coordinate of fluid points
void render_liquid(float *points, float diameter_pixels, int num_points, liquid_t *state)
{
    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    // Orphan current buffer
    glBufferData(GL_ARRAY_BUFFER, 2*num_points*sizeof(GLfloat), NULL, GL_STREAM_DRAW);

    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 2*num_points*sizeof(GLfloat), points, GL_STREAM_DRAW);

    // Unbind buffer
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Hack for reduced texture size
    float hack_diameter = diameter_pixels/(float)state->reduction;
    draw_liquid(state, hack_diameter, num_points);
}

void create_liquid_buffers(liquid_t *state)
{
    // VAO is REQUIRED for OpenGL 3+ when using VBO I believe
    #ifndef RASPI
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
    glGenBuffers(1, &state->tex_vbo);
    // Generate element buffer
    glGenBuffers(1, &state->tex_ebo);

    // Create frame buffer object for render to textures
    glGenFramebuffers(1, &state->frame_buffer_two);
    glBindFramebuffer(GL_FRAMEBUFFER, state->frame_buffer_two);

    // Create texture buffer
    glGenTextures(1, &state->tex_uniform);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);

    // http://fumufumu.q-games.com/gdc2010/shooterGDC.pdf
    // Render to low-res texture and upsample
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, state->screen_width/state->reduction, state->screen_height/state->reduction, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Attach image to framebuffer
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, state->tex_uniform, 0);

    #ifdef RASPI
    // Create frame buffer object for render to textures
    glGenFramebuffers(1, &state->frame_buffer_three);
    glBindFramebuffer(GL_FRAMEBUFFER, state->frame_buffer_three);
    #endif

    // Blur step requires additional texture to write into
    // Create texture buffer
    glGenTextures(1, &state->blur_horz_tex_uniform);
    glBindTexture(GL_TEXTURE_2D, state->blur_horz_tex_uniform);

    // http://fumufumu.q-games.com/gdc2010/shooterGDC.pdf
    // Render to low-res texture and upsample
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, state->screen_width/state->reduction, state->screen_height/state->reduction, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    #ifdef RASPI
    // Attach image to framebuffer
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, state->blur_horz_tex_uniform, 0);
    #else
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, state->blur_horz_tex_uniform, 0);
    #endif

    // Reset frame buffer and texture
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void create_texture_verticies(liquid_t *state)
{

    // Vertices: Pos(x,y) Tex(x,y)
    float vertices[] = {
        -1.0f, -1.0f, 0.0f, 1.0f,  // Bottom left
         1.0f, -1.0f, 1.0f, 1.0f, // Bottom right
         1.0f,  1.0f, 1.0f, 0.0f, // Top right
        -1.0f,  1.0f, 0.0f, 0.0f // Top left
    };

   // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->tex_vbo);
    // Fill buffer
    glBufferData(GL_ARRAY_BUFFER, 3*4*4*sizeof(GLfloat), vertices, GL_STATIC_DRAW);

    // Elements
    GLubyte elements[] = {
        2, 3, 0,
        0, 1, 2
    };

    // Set buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state->tex_ebo);
    // Fill buffer
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 2*3*sizeof(GLubyte), elements, GL_STATIC_DRAW);

}

void create_liquid_shaders(liquid_t *state)
{
    // Compile metaball vertex shader
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
      compile_shader(vertexShader, "SPH/shaders/liquid_particle_es.vert");
    #else
      compile_shader(vertexShader, "shaders/liquid_particle.vert");
    #endif

    // Compile metaball frag shader
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
      compile_shader(fragmentShader, "SPH/shaders/liquid_particle_es.frag");
    #else
      compile_shader(fragmentShader, "shaders/liquid_particle.frag");
    #endif

    // Create metaball shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertexShader);
    glAttachShader(state->program, fragmentShader); 

    // Link and use program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Compile texture vertex shader
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
      compile_shader(vertexShader, "SPH/shaders/render_liquid_texture_es.vert");
    #else
      compile_shader(vertexShader, "shaders/render_liquid_texture.vert");
    #endif

    // Compile texture frag shader
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
      compile_shader(fragmentShader, "SPH/shaders/render_liquid_texture_es.frag");
    #else
      compile_shader(fragmentShader, "shaders/render_liquid_texture.frag");
    #endif

    // Create texture shader program
    state->tex_program = glCreateProgram();
    glAttachShader(state->tex_program, vertexShader);
    glAttachShader(state->tex_program, fragmentShader);

    // Link and use program
    glLinkProgram(state->tex_program);
    show_program_log(state->tex_program);

    // Compile vert blur vertex shader
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
      compile_shader(vertexShader, "SPH/shaders/vert_blur_es.vert");
    #else
      compile_shader(vertexShader, "shaders/vert_blur.vert");
    #endif

    // Compile blur fragment shader(shared between horz vert blur vertex shaders)
    fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef RASPI
      compile_shader(fragmentShader, "SPH/shaders/blur_es.frag");
    #else
      compile_shader(fragmentShader, "shaders/blur.frag");
    #endif

    // Create vert blur shader program
    state->vert_blur_program = glCreateProgram();
    glAttachShader(state->vert_blur_program, vertexShader);
    glAttachShader(state->vert_blur_program, fragmentShader);

    // Link and use program
    glLinkProgram(state->vert_blur_program);
    show_program_log(state->vert_blur_program);    

    // Compile horz blur vertex shader
    vertexShader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef RASPI
      compile_shader(vertexShader, "SPH/shaders/horz_blur_es.vert");
    #else
      compile_shader(vertexShader, "shaders/horz_blur.vert");
    #endif

    // Create horz blur shader program
    state->horz_blur_program = glCreateProgram();
    glAttachShader(state->horz_blur_program, vertexShader);
    glAttachShader(state->horz_blur_program, fragmentShader);
    
    // Link and use program
    glLinkProgram(state->horz_blur_program);
    show_program_log(state->horz_blur_program);

    // Get position location
    state->position_location = glGetAttribLocation(state->program, "position");
    // Get pixel diameter location
    state->diameter_pixels_location = glGetUniformLocation(state->program, "diameter_pixels");

    // Get position location
    state->tex_position_location = glGetAttribLocation(state->tex_program, "position");
    // Get tex_coord location
    state->tex_coord_location = glGetAttribLocation(state->tex_program, "tex_coord");
    // Get tex uniform location
    state->tex_location = glGetUniformLocation(state->tex_program, "tex");

    // Get position location
    state->vert_blur_position_location = glGetAttribLocation(state->vert_blur_program, "position");
    // Get tex_coord location
    state->vert_blur_tex_coord_location = glGetAttribLocation(state->vert_blur_program, "tex_coord");
    // Get tex uniform location
    state->vert_blur_tex_location = glGetUniformLocation(state->vert_blur_program, "tex");

    // Get position location
    state->horz_blur_position_location = glGetAttribLocation(state->horz_blur_program, "position");
    // Get tex_coord location
    state->horz_blur_tex_coord_location = glGetAttribLocation(state->horz_blur_program, "tex_coord");
    // Get tex uniform location
    state->horz_blur_tex_location = glGetUniformLocation(state->horz_blur_program, "tex");

    // Enable point size to be specified in the shader
    #ifndef RASPI
    glEnable(GL_PROGRAM_POINT_SIZE);
    #endif

//   GLfloat fSizes[2];
//   glGetFloatv(GL_POINT_SIZE_RANGE,fSizes);
//   printf("min: %f, max: %f\n", fSizes[0], fSizes[1]);
}

void draw_liquid(liquid_t *state, float diameter_pixels, int num_points)
{
    //////
    // First phase - draw gaussian balls at particle position
    /////

    // Bind circle shader program
    glUseProgram(state->program);

    // Set pixel diameter uniform
    glUniform1f(state->diameter_pixels_location, (GLfloat)diameter_pixels);

    // Set buffer
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);

    glVertexAttribPointer(state->position_location, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GL_FLOAT), 0);
    glEnableVertexAttribArray(state->position_location);

    // Blend is required to show cleared color when the frag shader draws transparent pixels
    // Additive blending used for initial alpha draw
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE);

    // Bind frame buffer for render to texture
    glBindFramebuffer(GL_FRAMEBUFFER, state->frame_buffer_two);

    // Set viewport for low resolution texture
    glViewport(0,0,state->screen_width/state->reduction, state->screen_height/state->reduction);

    #ifndef RASPI
    // Set color attachment to draw into
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    #endif

    // Set background color
    glClearColor(0.0, 0.0, 0.0, 0.0);
    // Clear background
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw to color attachment 0 texture
    glDrawArrays(GL_POINTS, 0, num_points);

    //////
    // Second phase - horizontal blur
    /////
    // Bind vertical blur shader program
    glUseProgram(state->vert_blur_program);

    // Setup buffers
    size_t vert_size = 4*sizeof(GL_FLOAT);
    glBindBuffer(GL_ARRAY_BUFFER, state->tex_vbo);
    glVertexAttribPointer(state->vert_blur_position_location, 2, GL_FLOAT, GL_FALSE, vert_size, 0);
    glEnableVertexAttribArray(state->vert_blur_position_location);
    glVertexAttribPointer(state->vert_blur_tex_coord_location, 2, GL_FLOAT, GL_FALSE, vert_size,(void*)(2*sizeof(GL_FLOAT)));
    glEnableVertexAttribArray(state->vert_blur_tex_coord_location);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, state->tex_ebo);

    // Viewport already set for low rez texture

    // Do not need blending when drawing into 0,0,0,0 buffer
    glDisable(GL_BLEND);

    // Setup texture to read from
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->vert_blur_tex_location, 0);

    #ifdef RASPI
        // Bind frame buffer for render to texture
        glBindFramebuffer(GL_FRAMEBUFFER, state->frame_buffer_three);
    #else
        // Set color attachment to draw into
        glDrawBuffer(GL_COLOR_ATTACHMENT1);
    #endif

    // Set background color
    glClearColor(0.0, 0.0, 0.0, 0.0);
    // Clear background
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw to color attachment 1 texture
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

    //////
    // Third phase - vertical blur into original texture
    /////

    // Bind horizontal blur shader program
    glUseProgram(state->horz_blur_program);
    
    // Buffers already setup
    // Viewport already setup
    
    // Setup texture to read from
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->blur_horz_tex_uniform);
    glUniform1i(state->horz_blur_tex_location, 0);    


    #ifdef RASPI
    // Bind frame buffer for render to texture
    glBindFramebuffer(GL_FRAMEBUFFER, state->frame_buffer_two);
    #else
    // Set color attachment to draw into
    // Draw back into original attachment
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    #endif

    // Set background color
    glClearColor(0.0, 0.0, 0.0, 0.0);
    // Clear background
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw to color attachment 0 texture
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);

    //////
    // Fourth phase - Draw fluid image based on blured up sampled texture
    /////

    // Bind default frame buffer
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // Bind texture program
    glUseProgram(state->tex_program);

    // Buffers already setup

    // Set viewport back to "normal"
    glViewport(0,0,state->screen_width, state->screen_height);

    // Setup texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->tex_location, 0);

    // Enable "standard" blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Draw texture to screen
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);
}
