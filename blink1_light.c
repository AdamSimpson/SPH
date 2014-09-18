#include "blink1_light.h"
#include "blink1-lib.h"

// Set light to white
void rgb_light_white(rgb_light_t *state)
{
    blink1_setRGB( state->dev, 100, 100, 100);
}

// Reset color to state stored color
void rgb_light_reset(rgb_light_t *state)
{
    blink1_setRGB( state->dev, state->color[0], state->color[1], state->color[2] );
}

void init_rgb_light(rgb_light_t *state, uint8_t r, uint8_t g, uint8_t b) {
    state->color[0] = r;
    state->color[1] = g;
    state->color[2] = b;

    state->dev =  blink1_open();
    blink1_setRGB( state->dev, state->color[0], state->color[1], state->color[2] );
}

void shutdown_rgb_light(rgb_light_t *state) {
    blink1_setRGB( state->dev, 0, 0, 0);
    blink1_close(state->dev);
}
