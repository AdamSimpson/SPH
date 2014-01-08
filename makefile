CC=mpicc
CFLAGS= -mfloat-abi=hard -mfpu=vfp -O3 -lm

all:
	$(CC) $(CFLAGS) geometry.c hash.c fileio.c communication.c fluid.c -o sph.out
clean:
	rm -f ./*.o
