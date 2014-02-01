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
    font_state->screen_width = screen_width;
    font_state->screen_height = screen_height;

    // Create GL stash for 512x512 texture, our coordinate system has zero at bottom-left.
    font_state->fs = gl_shader_fonsCreate(512, 512, screen_width, screen_height, FONS_ZERO_BOTTOMLEFT);

    // Add font to stash.
    font_state->font_normal = fonsAddFont(font_state->fs, "sans", "DroidSerif-Regular.ttf");
    fonsSetFont(font_state->fs, font_state->font_normal);
    fonsSetSize(font_state->fs, 24.0f);
}

void render_fps(FONT_T *font_state, double fps)
{
    float dx, dy, lh;

    // Set font color
    unsigned int color = glfonsRGBA(255,255,100,255);
    fonsSetColor(font_state->fs, color);
    fonsVertMetrics(font_state->fs, NULL, NULL, &lh);

    // dx,dy,lh is internally handled by fontstash in pixels...arg
    dx = font_state->screen_width/2.0 - 5.0*lh;
    dy = font_state->screen_height/2.0 - lh;

    // Buffer to create strings in
    char buffer[64];

    // Gravity
    sprintf( buffer, "FPS: %.2f", fps);
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL); 
}

void render_parameters(FONT_T *font_state, double gravity, double viscocity, double density, double pressure, double elasticity)
{
    // Render some text
    float dx, dy, lh;
    unsigned int color = glfonsRGBA(255,255,100,255);

    // Get font height
    fonsVertMetrics(font_state->fs, NULL, NULL, &lh);
    fonsSetColor(font_state->fs, color);

    // dx,dy,lh is internally handled by fontstash in pixels...arg
    dx = -font_state->screen_width/2.0 + 10.0f;
    dy = font_state->screen_height/2.0 - lh;

    // Buffer to create strings in
    char buffer[100];

    // Gravity
    sprintf( buffer, "Gravity: %.1f", gravity);
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);

    // Viscocity
    sprintf( buffer, "Viscocity: %.1f", viscocity);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);

    // Density
    sprintf( buffer, "Density: %.1f", density);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);

    // Pressure
    sprintf( buffer, "Pressure: %.1f", pressure);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);

    // Elasticity
    sprintf( buffer, "Elasticity: %.1f", elasticity);
    dy -= lh;
    fonsDrawText(font_state->fs, dx, dy, buffer, NULL);

}

void remove_font(FONT_T *font_state)
{
    gl_shader_fonsDelete(font_state->fs);

}
