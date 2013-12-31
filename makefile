CC=cc
CFLAGS=-O3 -lm

all:
	$(CC) $(CFLAGS) geometry.c hash.c fileio.c fluid.c -o sph.out
