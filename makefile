CC=cc
CFLAGS=-g -lm

all:
	$(CC) $(CFLAGS) hash.c fileio.c fluid.c -o sph.out
