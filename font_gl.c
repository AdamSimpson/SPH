#include <stdio.h>
#include <stdlib.h>

#define FONTSTASH_IMPLEMENTATION
#include "fontstash.h"

#include "fluid.h"
#include "font_gl.h"

#define GL_SHADER_FONTSTASH_IMPLEMENTATION
#include "gl_shader_fontstash.h"

void init_font(FONT_T *font_state, int screen_width, int screen_height)
{
    // Create GL stash for 512x512 texture, our coordinate system has zero at bottom-left.
    font_state->fs = gl_shader_fonsCreate(512, 512, screen_width, screen_height, FONS_ZERO_BOTTOMLEFT);

    // Add font to stash.
    font_state->font_normal = fonsAddFont(font_state->fs, "sans", "DroidSerif-Regular.ttf");
}

void render_font(FONT_T *font_state)
{
    // Render some text
    float dx, dy, lh;
    unsigned int color = glfonsRGBA(255,255,100,255);

    // Get font height
    fonsVertMetrics(font_state->fs, NULL, NULL, &lh);

    fonsSetFont(font_state->fs, font_state->font_normal);
    fonsSetSize(font_state->fs, 64.0f);
    fonsSetColor(font_state->fs, color);

    dx = -400.0;
    dy = 0.0;
    dx = fonsDrawText(font_state->fs, dx, dy,"Gravity: ", NULL);
    fonsDrawText(font_state->fs, dx, dy, "+1", NULL);
}

void remove_font(FONT_T *font_state)
{
    gl_shader_fonsDelete(font_state->fs);

}
