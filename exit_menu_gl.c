#include <stdio.h>
#include <unistd.h>
#include "exit_menu_gl.h"
#include "renderer.h"
#include "image_gl.h"
#include "cursor_gl.h"

void init_exit_menu(exit_menu_t *state, gl_t *gl_state)
{
    // Initialize program launcher images OpenGL state
    int image_width = 500;
    int image_height = 312;

    float half_height = (image_height/(float)gl_state->screen_height);
    float half_width = (image_width/(float)gl_state->screen_width);

    float dx = (2.0f-(3.0*2.0f*half_width))/4.0f;

    float lower_left_y = -half_height;
    float lower_left_x = -1.0f + dx;
    state->mandelbrot_state = malloc(sizeof(image_t));

    #ifdef RASPI
    init_image(state->mandelbrot_state,
               gl_state,
               "SPH/mandelbrot.png",
               "SPH/mandelbrot-selected.png",
               lower_left_x, lower_left_y,
               image_width, image_height);
    #else
    init_image(state->mandelbrot_state,
               gl_state,
               "mandelbrot.png",
               "mandelbrot-selected.png",
               lower_left_x, lower_left_y,
               image_width, image_height);
    #endif

    lower_left_x += 2.0f*half_width + dx;
    state->terminal_state = malloc(sizeof(image_t));
    #ifdef RASPI
    init_image(state->terminal_state,
               gl_state,
               "SPH/terminal.png",
               "SPH/terminal-selected.png",
               lower_left_x, lower_left_y,
               image_width, image_height);
    #else
    init_image(state->terminal_state,
               gl_state,
               "terminal.png",
               "terminal-selected.png",
               lower_left_x, lower_left_y,
               image_width, image_height);
    #endif

    lower_left_x += 2.0f*half_width + dx;
    state->sph_state = malloc(sizeof(image_t));
    #ifdef RASPI
    init_image(state->sph_state,
               gl_state,
               "SPH/sph.png",
               "SPH/sph-selected.png",
               lower_left_x, lower_left_y,
               image_width, image_height);
    #else
    init_image(state->sph_state,
               gl_state,
               "sph.png",
               "sph-selected.png",
               lower_left_x, lower_left_y,
               image_width, image_height);
    #endif

    // Initialize cursor
    state->cursor_state = malloc(sizeof(cursor_t));
    int cursor_width = 30;
    int cursor_height = 40;
    #ifdef RASPI
    init_cursor(state->cursor_state, gl_state, "SPH/cursor.png", cursor_width, cursor_height);
    #else
    init_cursor(state->cursor_state, gl_state, "cursor.png", cursor_width, cursor_height); 
    #endif
}

void exit_exit_menu(exit_menu_t *state)
{
    free(state->mandelbrot_state);
    free(state->sph_state);
    free(state->terminal_state);
    free(state->cursor_state);
}

void render_exit_menu(exit_menu_t *state, float cursor_x, float cursor_y)
{
    // Update center of cursor
    set_cursor_position(state->cursor_state, cursor_x, cursor_y);

    // Check if anything selected
    check_cursor_in_image(state->cursor_state, state->mandelbrot_state);
    check_cursor_in_image(state->cursor_state, state->sph_state);
    check_cursor_in_image(state->cursor_state, state->terminal_state);

    // Draw mandelbrot image
    draw_image(state->mandelbrot_state);

    // Draw terminal image
    draw_image(state->terminal_state);

    // Draw SPH image
    draw_image(state->sph_state);

    // Draw cursor
    draw_cursor(state->cursor_state);
}

void update_cursor(exit_menu_t *state, float x, float y)
{
    // Update cursor state
    state->cursor_state->center_x = x;
    state->cursor_state->center_y = y;
}

// Check if cursor is inside of image, if so set selected attribute
void check_cursor_in_image(cursor_t *cursor_state, image_t *image_state)
{
    float x, y;

    // image dimensions in gl screen coordinates
    float half_width = (cursor_state->cursor_width/(float)cursor_state->gl_state->screen_width);
    float half_height = (cursor_state->cursor_height/(float)cursor_state->gl_state->screen_height);

    // Upper left corner
    x = cursor_state->center_x - half_width;
    y = cursor_state->center_y + half_height;
    if( inside_image(image_state, x, y) ){
        image_state->selected = true;
        return;
    }
    // Lower left corner
    x = cursor_state->center_x - half_width;
    y = cursor_state->center_y - half_height;
    if( inside_image(image_state, x, y) ){
        image_state->selected = true;
        return;
    }
    // Upper right corner
    x = cursor_state->center_x + half_width;
    y = cursor_state->center_y + half_height;
    if( inside_image(image_state, x, y) ){
        image_state->selected = true;
        return;
    }
    // Lower right corner
    x = cursor_state->center_x + half_width;
    y = cursor_state->center_y - half_height;
    if( inside_image(image_state, x, y) ){
        image_state->selected = true;
        return;
    }

    image_state->selected = false;
}

