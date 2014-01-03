CC=cc
CFLAGS= -O3 -lm

all:
	$(CC) $(CFLAGS) geometry.c hash.c fileio.c communication.c fluid.c -o sph.out
	cp sph.out /tmp/work/${USER}
