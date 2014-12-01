#include "rgb_light.h"

#define ARRAY_SIZE(a) (sizeof(a) / sizeof((a)[0]))

const uint8_t gamma[] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2,
  2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5,
  5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 9, 9, 9, 10,
  10, 10, 11, 11, 11, 12, 12, 13, 13, 13, 14, 14, 15, 15, 16, 16,
  17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 24, 24, 25,
  25, 26, 27, 27, 28, 29, 29, 30, 31, 32, 32, 33, 34, 35, 35, 36,
  37, 38, 39, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 50,
  51, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 66, 67, 68,
  69, 70, 72, 73, 74, 75, 77, 78, 79, 81, 82, 83, 85, 86, 87, 89,
  90, 92, 93, 95, 96, 98, 99,101,102,104,105,107,109,110,112,114,
  115,117,119,120,122,124,126,127,129,131,133,135,137,138,140,142,
  144,146,148,150,152,154,156,158,160,162,164,167,169,171,173,175,
  177,180,182,184,186,189,191,193,196,198,200,203,205,208,210,213,
  215,218,220,223,225,228,231,233,236,239,241,244,247,249,252,255 };

void pabort(const char *s)
{
    perror(s);
    abort();
}

// Set light color off
void rgb_light_off(rgb_light_t *state)
{
    transfer(state, 0, 0, 0);
}

// Set light to white
void rgb_light_white(rgb_light_t *state)
{
    transfer(state, 4, 4, 4);
}

// Reset color to state stored color
void rgb_light_reset(rgb_light_t *state)
{
    transfer(state, state->color[0], state->color[1], state->color[2]);
}

void transfer(rgb_light_t *state, uint8_t r, uint8_t g, uint8_t b)
{
    int ret;
    uint8_t tx[] = {gamma[r], gamma[g], gamma[b]};
    uint8_t rx[ARRAY_SIZE(tx)] = {0, };
    struct spi_ioc_transfer tr = {
      .tx_buf = (unsigned long)tx,
      .rx_buf = (unsigned long)rx,
      .len = ARRAY_SIZE(tx),
      .delay_usecs = state->delay,
      .speed_hz = state->speed,
      .bits_per_word = state->bits,
    };

    ret = ioctl(state->fd, SPI_IOC_MESSAGE(1), &tr);
    if (ret < 1)
      pabort("can't send spi message");
}

void init_rgb_light(rgb_light_t *state, uint8_t r, uint8_t g, uint8_t b) {
    int ret = 0;

    state->mode = 0;
    state->bits = 8;
    state->speed = 2000000;
    state->delay = 0;
    state->fd = open("/dev/spidev0.0", O_RDWR);
    state->color[0] = r;
    state->color[1] = g;
    state->color[2] = b;
    if (state->fd < 0)
        pabort("can't open device");

    /*
     * spi mode
     */
    ret = ioctl(state->fd, SPI_IOC_WR_MODE, &state->mode);
    if (ret == -1)
        pabort("can't set spi mode");

    ret = ioctl(state->fd, SPI_IOC_RD_MODE, &state->mode);
    if (ret == -1)
        pabort("can't get spi mode");

    /*
     * bits per word
     */
    ret = ioctl(state->fd, SPI_IOC_WR_BITS_PER_WORD, &state->bits);
    if (ret == -1)
        pabort("can't set bits per word");

    ret = ioctl(state->fd, SPI_IOC_RD_BITS_PER_WORD, &state->bits);
    if (ret == -1)
        pabort("can't get bits per word");

    /*
     * max speed hz
     */
    ret = ioctl(state->fd, SPI_IOC_WR_MAX_SPEED_HZ, &state->speed);
      if (ret == -1)
          pabort("can't set max speed hz");

    ret = ioctl(state->fd, SPI_IOC_RD_MAX_SPEED_HZ, &state->speed);
    if (ret == -1)
        pabort("can't get max speed hz");

    /* printf("spi mode: %d\n", mode); */
    /* printf("bits per word: %d\n", bits); */
    /* printf("max speed: %d Hz (%d KHz)\n", speed, speed/1000); */

    transfer(state, state->color[0], state->color[1], state->color[2]);
}

void shutdown_rgb_light(rgb_light_t *state) {
  rgb_light_off(state);
  close(state->fd);
}
