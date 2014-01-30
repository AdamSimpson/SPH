#include <stdio.h>
#include <stdlib.h>

#define FONTSTASH_IMPLEMENTATION
#include "fontstash.h"

#include "font_gl.h"

#define GL_SHADER_FONTSTASH_IMPLEMENTATION
#include "gl_shader_fontstash.h"

void init_font(FONT_T *font_state)
{
    // Create GL stash for 512x512 texture, our coordinate system has zero at top-left.
    font_state->fs = gl_shader_fonsCreate(512, 512, FONS_ZERO_TOPLEFT);

    // Add font to stash.
    font_state->font_normal = fonsAddFont(font_state->fs, "sans", "DroidSerif-Regular.ttf");
}

void render_font(FONT_T *font_state)
{

    // Render some text
    float dx = 10, dy = 10;
    unsigned int white = glfonsRGBA(255,255,255,255);
    unsigned int brown = glfonsRGBA(192,128,0,128);

    fonsSetFont(font_state->fs, font_state->font_normal);
    fonsSetSize(font_state->fs, 124.0f);
    fonsSetColor(font_state->fs, white);
    fonsDrawText(font_state->fs, dx,dy,"SPH", NULL);
}

void remove_font(FONT_T *font_state)
{
    gl_shader_fonsDelete(font_state->fs);

}
