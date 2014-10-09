#ifndef RGB_LIGHT_H
#define RGB_LIGHT_H

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "blink1-lib.h"

typedef struct rgb_light_t {
    blink1_device* dev; 
    uint8_t color[3];
} rgb_light_t;

void rgb_light_white(rgb_light_t *state);
void rgb_light_reset(rgb_light_t *state);
void init_rgb_light(rgb_light_t *state, uint8_t r, uint8_t g, uint8_t b);
void shutdown_rgb_light(rgb_light_t *state);

#endif
