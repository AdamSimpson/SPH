CC=mpicc
CFLAGS=-lm -g -O3

all:
	$(CC) $(CFLAGS) geometry.c hash.c fileio.c communication.c fluid.c -o sph.out
clang:
	$(CC) -compile_info -link_info -Wl,-flat_namespace -lm -g geometry.c hash.c fileio.c communication.c fluid.c -o sph.out -I/usr/local/Cellar/mpich2/3.1.3_1/include -L/usr/local/Cellar/mpich2/3.1.3_1/lib -lmpi -lpmpi

clean:
	rm -f ./*.o
