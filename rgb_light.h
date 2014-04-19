#ifndef RGB_LIGHT_H
#define RGB_LIGHT_H

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <linux/types.h>
#include <linux/spi/spidev.h>
#include <signal.h>

typedef struct rgb_light_t {
    uint8_t mode;
    uint8_t bits;
    uint32_t speed;
    uint16_t delay;
    int fd;
    uint8_t color[3];
} rgb_light_t;

void pabort(const char *s);
void transfer(rgb_light_t *state, uint8_t r, uint8_t g, uint8_t b);
void rgb_light_off(rgb_light_t *state);
void rgb_light_white(rgb_light_t *state);
void rgb_light_reset(rgb_light_t *state);
void init_rgb_light(rgb_light_t *state, uint8_t r, uint8_t g, uint8_t b);
void rgb_light_off(rgb_light_t *state);
void shutdown_rgb_light(rgb_light_t *state);

#endif
