CC=mpicc
CFLAGS=-O3 -lm -g

all:
	$(CC) $(CFLAGS) geometry.c hash.c fileio.c communication.c fluid.c -o sph.out
clean:
	rm -f ./*.o
