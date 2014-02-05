#include <stdio.h>
#include <stdlib.h>
#include "fluid.h"
#include "renderer.h"
#include "font_gl.h"

#include "ogl_utils.h"

#include <ft2build.h>
#include FT_FREETYPE_H

#define MAX_WIDTH 2048 // Maximum texture width on pi

void create_font_program(FONT_T *state)
{
    // Compile vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    #ifdef GLES
        compile_shader(vertex_shader, "SPH/font_es.vert");
    #else
        compile_shader(vertex_shader, "font.vert");
    #endif

    // Compile fragment shader
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    #ifdef GLES
        compile_shader(frag_shader, "SPH/font_es.frag");
    #else
        compile_shader(frag_shader, "font.frag");
    #endif

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertex_shader);
    glAttachShader(state->program, frag_shader);

    // Link  program
    glLinkProgram(state->program);
    show_program_log(state->program);


    // Get coord attribute location
    state->coord_location = glGetAttribLocation(state->program, "coord");
    // Get tex uniform location
    state->tex_location = glGetUniformLocation(state->program, "tex");
    // Get color uniform location
    state->color_location = glGetUniformLocation(state->program, "color");

    // Enable blend
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void create_font_buffers(FONT_T *state)
{
    // VAO is required for OpenGL 3+ when using VBO I believe
    #ifndef GLES
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    #endif

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

void create_font_atlas(FONT_T *state)
{
    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &state->tex_uniform);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->tex_uniform, 0);

    // Set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Set single byte alignment
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    // Get atlas dimensions
    FT_GlyphSlot g = state->face->glyph;
    int w = 0; // full texture width
    int h = 0; // full texture height
    int row_w = 0; // current row width
    int row_h = 0; // current row height

    int i;
    for(i=32; i<128; i++) {
        if(FT_Load_Char(state->face, i, FT_LOAD_RENDER)) {
            printf("Loading Character %d failed\n", i);
            exit(EXIT_FAILURE);
        }

        // If the width will be over max texture width
        // Go to next row
        if(row_w + g->bitmap.width+1 >= MAX_WIDTH) {
            w = max(w, row_w);
            h += row_h;
            row_w = 0;
            row_h = 0;
        }
        row_w += g->bitmap.width + 1;
        row_h = max(row_h, g->bitmap.rows);
    }
    
    // final texture dimensions
    w = max(row_w, w);
    h += row_h;


    state->atlas_width = w;
    state->atlas_height = h;

    // Allocate texture
    #ifdef GLES
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, w, h, 0, GL_ALPHA, GL_UNSIGNED_BYTE, 0);
    #else
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_UNSIGNED_BYTE, 0);
    #endif

    // Fill texture with glyph bitmaps and cache placements
    CHAR_INFO *char_info = state->char_info;
    int offset_x = 0;
    int offset_y = 0;
    row_h = 0;

    for(i=32; i<128; i++) {
        if(FT_Load_Char(state->face, i, FT_LOAD_RENDER)) {
            printf("Loading Character %d failed\n", i);
            exit(EXIT_FAILURE);
        }

        // Set correct row
        if(offset_x + g->bitmap.width + 1 >= MAX_WIDTH) {
            offset_y += row_h;
            row_h = 0;
            offset_x = 0;
        }

        // fill texture with glyph
        glTexSubImage2D(GL_TEXTURE_2D, 0, offset_x, offset_y, g->bitmap.width, g->bitmap.rows, GL_ALPHA, GL_UNSIGNED_BYTE, g->bitmap.buffer);

        // Cache values
        char_info[i].ax = g->advance.x >> 6;
        char_info[i].ay = g->advance.y >> 6;
        char_info[i].bw = g->bitmap.width;
        char_info[i].bh = g->bitmap.rows;
        char_info[i].bl = g->bitmap_left;
        char_info[i].bt = g->bitmap_top;
        char_info[i].tx = offset_x/(float)w;
        char_info[i].ty = offset_y/(float)h;

        // Update current position
        row_h = max(row_h, g->bitmap.rows);
        offset_x += g->bitmap.width + 1;
    }
} 

void render_text(FONT_T *state, const char *text, float x, float y, float sx, float sy)
{
    // Setup environment
    glUseProgram(state->program);
    glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
    glVertexAttribPointer(state->coord_location, 4, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(state->coord_location);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
    glUniform1i(state->tex_location, 0);

    // Set font color
    GLfloat black[4] = {1, 1, 1, 1};
    glUniform4fv(state->color_location, 1, black);

    struct point {
        GLfloat x;
        GLfloat y;
        GLfloat s;
        GLfloat t;
    } coords[6*strlen(text)];

    int n = 0;

    CHAR_INFO *c = state->char_info;
    
    char *p;
    for(p = text; *p; p++) {
        float x2 = x + c[*p].bl * sx;
        float y2 = -y - c[*p].bt * sy;
        float w = c[*p].bw * sx;
        float h = c[*p].bh * sy;

        // Advance cursor to start of next char
        x += c[*p].ax * sx;
        y += c[*p].ay * sy;

        // Skip 0 pixel glyphs
        if(!w || !h)
            continue;

        coords[n++] = (struct point){x2, -y2, c[*p].tx, c[*p].ty};
        coords[n++] = (struct point){x2+w, -y2, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty};
        coords[n++] = (struct point){x2, -y2-h, c[*p].tx, c[*p].ty + c[*p].bh/state->atlas_height};
        coords[n++] = (struct point){x2+w, -y2, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty};
        coords[n++] = (struct point){x2, -y2-h, c[*p].tx, c[*p].ty + c[*p].bh/state->atlas_height};
        coords[n++] = (struct point){x2 + w, -y2-h, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty + c[*p].bh/state->atlas_height};

        glBufferData(GL_ARRAY_BUFFER, sizeof(coords), coords, GL_DYNAMIC_DRAW);
        glDrawArrays(GL_TRIANGLES, 0, n);
    }
}

void init_font(FONT_T *state, int screen_width, int screen_height)
{
    // Initialize FreeType library
    if(FT_Init_FreeType(&state->ft)) {
        printf("Error initializing FreeType library\n");
        exit(EXIT_FAILURE);
    }

    state->screen_width = screen_width;
    state->screen_height = screen_height;

    // Load font face
    if(FT_New_Face(state->ft, "SPH/DroidSerif-Regular.ttf", 0, &state->face)) {
        printf("Error loading font face\n");
        exit(EXIT_FAILURE);
    }

    // Set pixel size
    FT_Set_Pixel_Sizes(state->face, 0, 24);

    // Setup OpenGL
    create_font_program(state);
    create_font_buffers(state);
    create_font_atlas(state);
}

void render_fps(FONT_T *state, double fps)
{
    float dx, dy, lh;

    // Font start
    float sx = 2.0 / state->screen_width;
    float sy = 2.0 / state->screen_height;

    // Buffer to create strings in
    char buffer[64];
    sprintf( buffer, "FPS: %.2f", fps);
   
    // Render text
    render_text(state, buffer, -1 + 8 * sx, 1 - 50 * sy, sx, sy);

}

void render_parameters(FONT_T *font_state, parameters selected_param, double gravity, double viscosity, double density, double pressure, double elasticity)
{
/*
    // Render some text
    float dx, dy, lh;
    unsigned int color = glfonsRGBA(255,255,100,255);

    // Get font height
    fonsVertMetrics(font_state->fs, NULL, NULL, &lh);

    // dx,dy,lh is internally handled by fontstash in pixels...arg
    dx = -font_state->screen_width/2.0 + 10.0f;
    dy = font_state->screen_height/2.0 - lh;

    // Buffer to create strings in
    char buffer[100];

    // Add blur to parameters
    //unblur selected parameter below
    fonsSetBlur(font_state->fs, 2.0f);

    // Gravity
    if(selected_param == GRAVITY)
        fonsSetBlur(font_state->fs, 0);
    sprintf( buffer, "Gravity: %.1f", gravity);
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);
    fonsSetBlur(font_state->fs, 2.0f);

    // Viscocity
    if(selected_param == VISCOSITY)
        fonsSetBlur(font_state->fs, 0);
    sprintf( buffer, "Viscosity: %.1f", viscosity);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);
    fonsSetBlur(font_state->fs, 2.0f);

    // Density
    if(selected_param == DENSITY)
        fonsSetBlur(font_state->fs, 0);
    sprintf( buffer, "Density: %.1f", density);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);
    fonsSetBlur(font_state->fs, 2.0f);

    // Pressure
    if(selected_param == PRESSURE)
        fonsSetBlur(font_state->fs, 0);
    sprintf( buffer, "Pressure: %.1f", pressure);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);
    fonsSetBlur(font_state->fs, 2.0f);

    // Elasticity
    if(selected_param == ELASTICITY)
        fonsSetBlur(font_state->fs, 0);
    sprintf( buffer, "Elasticity: %.1f", elasticity);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);


    // Unset blur
    fonsSetBlur(font_state->fs, 0);  
*/
}

void remove_font(FONT_T *font_state)
{

}
