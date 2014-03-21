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
    #ifdef GLES
        compile_shader(vertex_shader, "SPH/shaders/background_es.vert");
    #else
        compile_shader(vertex_shader, "shaders/background.vert");
    #endif

    // Compile fragment shader
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef GLES
        compile_shader(frag_shader, "SPH/background_es.frag");
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
    // Get tec_coord location
    state->tex_coord_location = glGetAttribLocation(state->program, "tex_coord");
    // Get tex uniform location
    state->tex_location = glGetUniformLocation(state->program, "tex");
}

void create_background_buffers(background_t *state)
{
    // VAO is required for OpenGL 3+ when using VBO I believe
    #ifndef GLES
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
    float vertices[] = {
         // Full screen vertices
        -1.0f,  1.0f, 0.0f, 0.0f, // Top left
         1.0f,  1.0f, 1.0f, 0.0f, // Top right
         1.0f, -1.0f, 1.0f, 1.0f, // Bottom right
	-1.0f, -1.0f, 0.0f, 1.0f  // Bottom left
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

void create_background_texture(background_t *state)
{
    // Read in background PNG
    unsigned error;
    unsigned char* image;
    unsigned width, height;

    error = lodepng_decode32_file(&image, &width, &height, "Titan_bgnd.png");
    if(error) printf("error %u: %s\n", error, lodepng_error_text(error));

    printf("Background image loaded: %d x %d pixels\n",width, height);

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
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    // Release image host memory
    free(image);
}

void init_background(background_t *state)
{
    // Create program
    create_backround_program(state);

    // Generate buffers 
    create_background_buffers(state);  

    // Set verticies
    create_background_vertices(state);

    // Create texture from image
    create_background_texture(state);
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
    glDisable(GL_BLEND);

    // Setup texture
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->tex_location, 0);

    // Draw background
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_BYTE, 0);
}
