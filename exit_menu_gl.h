#ifndef EXIT_H
#define EXIT_H

typedef struct EXIT_MENU_T exit_menu_t;

#include "cursor_gl.h"
#include "image_gl.h"

struct EXIT_MENU_T {
    image_t *mandelbrot_state;
    image_t *sph_state;
    image_t *terminal_state;
    cursor_t *cursor_state;
};

void exit_exit_menu(exit_menu_t *state);
void init_exit_menu(exit_menu_t *state, gl_t *gl_state);
void render_exit_menu(exit_menu_t *state, float cursor_x, float cursor_y);
void update_cursor(exit_menu_t *state, float x, float y);
void check_cursor_in_image(cursor_t *cursor_state, image_t *image_state);

#endif             
