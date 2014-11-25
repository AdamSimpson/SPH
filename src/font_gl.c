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
#include "font_gl.h"

#include "ogl_utils.h"

#include <ft2build.h>
#include FT_FREETYPE_H

#define MAX_WIDTH 2048 // Maximum texture width on pi

void create_font_program(font_t *state)
{
    // Compile vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    compile_shader(vertex_shader, "shaders/font.vert");

    // Compile fragment shader
    GLuint frag_shader = glCreateShader(GL_FRAGMENT_SHADER);
    compile_shader(frag_shader, "shaders/font.frag");

    // Create shader program
    state->program = glCreateProgram();
    glAttachShader(state->program, vertex_shader);
    glAttachShader(state->program, frag_shader);

    // Link  program
    glLinkProgram(state->program);
    show_program_log(state->program);

    // Get coord attribute location
    state->coord_location = glGetAttribLocation(state->program, "coord");
    // Get color attribute location
    state->color_location = glGetAttribLocation(state->program, "color");
    // Get tex uniform location
    state->tex_location = glGetUniformLocation(state->program, "tex");

    // Enable blend
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}

void create_font_buffers(font_t *state)
{
    // VAO is required for OpenGL 3+ when using VBO I believe
    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    // Generate vertex buffer
    glGenBuffers(1, &state->vbo);
}

// Create font atlas texture
void create_font_atlas(font_t *state)
{
    glUseProgram(state->program);
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
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, w, h, 0, GL_RED, GL_UNSIGNED_BYTE, 0);

    // Fill texture with glyph bitmaps and cache placements
    char_info_t *char_info = state->char_info;
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
	    glTexSubImage2D(GL_TEXTURE_2D, 0, offset_x, offset_y, g->bitmap.width, g->bitmap.rows, GL_RED, GL_UNSIGNED_BYTE, g->bitmap.buffer);	

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

// Render single string
// Assumes arrays and program already steup
void render_text(font_t *state, char *text, float x, float y, float sx, float sy)
{
	text_vert_t verts[6*strlen(text)];

	int n = 0;

	char_info_t *c = state->char_info;

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

		verts[n++] = (text_vert_t){x2, -y2, c[*p].tx, c[*p].ty, 1, 1, 1};
		verts[n++] = (text_vert_t){x2+w, -y2, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty, 1, 1, 1};
		verts[n++] = (text_vert_t){x2, -y2-h, c[*p].tx, c[*p].ty + c[*p].bh/state->atlas_height, 1, 1, 1};
		verts[n++] = (text_vert_t){x2+w, -y2, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty, 1, 1, 1};
		verts[n++] = (text_vert_t){x2, -y2-h, c[*p].tx, c[*p].ty + c[*p].bh/state->atlas_height, 1, 1, 1};
		verts[n++] = (text_vert_t){x2 + w, -y2-h, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty + c[*p].bh/state->atlas_height, 1, 1, 1};
	}

	// Orphan buffer
	glBufferData(GL_ARRAY_BUFFER, sizeof(verts), NULL, GL_STREAM_DRAW);

	// Buffer vertices
	glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STREAM_DRAW);

	// Draw text
	glDrawArrays(GL_TRIANGLES, 0, n);
}

// Add text coordinates to be rendered later
// This allows multiple strings to be rendered with a single buffer and draw
int add_text_coords(font_t *state, char *text, text_vert_t* verts, float *color, float x, float y, float sx, float sy)
{
	int n = 0;

	char_info_t *c = state->char_info;

	float r = color[0];
	float g = color[1];
	float b = color[2];

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

		verts[n++] = (text_vert_t){x2, -y2, c[*p].tx, c[*p].ty, r, g, b};
		verts[n++] = (text_vert_t){x2+w, -y2, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty, r, g, b};
		verts[n++] = (text_vert_t){x2, -y2-h, c[*p].tx, c[*p].ty + c[*p].bh/state->atlas_height, r, g, b};
		verts[n++] = (text_vert_t){x2+w, -y2, c[*p].tx + c[*p].bw/state->atlas_width, c[*p].ty, r, g, b};
		verts[n++] = (text_vert_t){x2, -y2-h, c[*p].tx, c[*p].ty + c[*p].bh/state->atlas_height, r, g, b};
		verts[n++] = (text_vert_t){x2 + w, -y2-h, c[*p].tx + c[*p].bw/state->atlas_width, 
			c[*p].ty + c[*p].bh/state->atlas_height, r, g, b};
	}

	return n;
}

void render_all_text(font_t *state, render_t *render_state, double fps)
{
	// Setup environment
	glUseProgram(state->program);
	glBindBuffer(GL_ARRAY_BUFFER, state->vbo);
	glVertexAttribPointer(state->coord_location, 4, GL_FLOAT, GL_FALSE, 7*sizeof(GLfloat), 0);
	glVertexAttribPointer(state->color_location, 3, GL_FLOAT, GL_FALSE, 7*sizeof(GLfloat), (void*)(4*sizeof(GLfloat)));
	glEnableVertexAttribArray(state->coord_location);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, state->tex_uniform);
	glUniform1i(state->tex_location, 0);

	// Blend is required to show cleared color when the frag shader draws transparent pixels
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// Buffer to hold individual string
	char buffer[100];

	// Text points
	// Max of 512 characters in total
	text_vert_t verts[6*512];

	// Font start
	float sx = 2.0f / state->screen_width;
	float sy = 2.0f / state->screen_height;

	parameters selected_param = render_state->selected_parameter;
	float gravity, viscosity, density, k, dq, smooth;

	gravity = render_state->master_params[0].g;
	viscosity = render_state->master_params[0].c;
	density = render_state->master_params[0].rest_density;
	k = render_state->master_params[0].k;
	dq = render_state->master_params[0].dq;
        smooth = render_state->master_params[0].smoothing_radius;

	float unselected_color[3] = {1.0f ,1.0f, 1.0f};
	float selected_color[3]   = {0.1f, 0.80f, 0.43f};  
	float *color;

	// Total number of font coordinates
	int n = 0;

	// frames per second
	sprintf( buffer, "FPS: %.0f", fps);
	n += add_text_coords(state, buffer, verts + n, unselected_color, 1.0f - 100.0f * sx, 1.0f - 50.0f * sy, sx, sy);

	// Gravity
	sprintf( buffer, "Gravity: %.1f", gravity);
	if(selected_param == GRAVITY)
		color = selected_color;
	else
		color = unselected_color;
	n += add_text_coords(state, buffer, verts + n, color, -1.0f + 8.0f * sx, 1.0f - 100.0f * sy, sx, sy);

        // Smoothing
        sprintf( buffer, "smooth: %.2f", smooth);
        if(selected_param == SMOOTH)
                color = selected_color;
        else
                color = unselected_color;
        n += add_text_coords(state, buffer, verts + n, color, -1.0f + 8.0f * sx, 1.0f - 150.0f * sy, sx, sy);

	// Density
	sprintf( buffer, "Density: %.2f", density);
	if(selected_param == DENSITY)
		color = selected_color;
	else
		color = unselected_color;
	n += add_text_coords(state, buffer, verts + n, color, -1.0f + 8.0f * sx, 1.0f - 200.0f * sy, sx, sy);

	// K
	sprintf( buffer, "K: %.2f", k);
	if(selected_param == K)
		color = selected_color;
	else
		color = unselected_color;
	n += add_text_coords(state, buffer, verts + n, color, -1.0f + 8.0f * sx, 1.0f - 250.0f * sy, sx, sy);

	// DQ
	sprintf( buffer, "dq: %.2f", dq);
	if(selected_param == DQ)
		color = selected_color;
	else
		color = unselected_color;
	n += add_text_coords(state, buffer, verts + n, color, -1.0f + 8.0f * sx, 1.0f - 300.0f * sy, sx, sy);

       // Viscocity
        sprintf( buffer, "Viscosity: %.2f", viscosity);
        if(selected_param == VISCOSITY)
                color = selected_color;
        else
                color = unselected_color;
        n += add_text_coords(state, buffer, verts + n, color, -1.0f + 8.0f * sx, 1.0f - 350.0f * sy, sx, sy);


	// Orphan buffer
	glBufferData(GL_ARRAY_BUFFER, n*sizeof(text_vert_t), NULL, GL_STREAM_DRAW);

	// Buffer vertices
	glBufferData(GL_ARRAY_BUFFER, n*sizeof(text_vert_t), verts, GL_STREAM_DRAW);

	// Draw text
	glDrawArrays(GL_TRIANGLES, 0, n);
}

// Setup freetype font
void init_font(font_t *state, int screen_width, int screen_height)
{
	// Initialize FreeType library
	if(FT_Init_FreeType(&state->ft)) {
		printf("Error initializing FreeType library\n");
		exit(EXIT_FAILURE);
	}

	state->screen_width = screen_width;
	state->screen_height = screen_height;

	// Load font face
	if(FT_New_Face(state->ft, "DroidSerif-Regular.ttf", 0, &state->face)) {
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

void remove_font(font_t *font_state)
{

}
