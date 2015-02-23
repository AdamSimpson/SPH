CC=mpicc
CFLAGS=-g -O3
LDLIBS=-lm

SRCS=geometry.c hash.c fileio.c communication.c fluid.c
OBJS=$(subst .c,.o, $(SRCS))

all: sph

sph: $(OBJS)
	$(CC) $(CFLAGS) $(LDLIBS) -o sph $(OBJS)

clang:
	$(CC) -compile_info -link_info -Wl,-flat_namespace -lm -g geometry.c hash.c fileio.c communication.c fluid.c -o sph.out -I/usr/local/Cellar/mpich2/3.1.3_1/include -L/usr/local/Cellar/mpich2/3.1.3_1/lib -lmpi -lpmpi

clean:
	$(RM) $(OBJS)
	$(RM) sph
