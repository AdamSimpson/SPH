CC=cc
CFLAGS=-g

all:
	$(CC) $(CFLAGS) geometry.c hash.c fileio.c fluid.c -o sph.out
