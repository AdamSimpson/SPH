CC=mpicc
CFLAGS=-g -O3
LDLIBS=-lm

SRC_DIR=../src
BIN_DIR=../bin
SRCS=geometry.c hash.c fileio.c communication.c fluid.c
OBJS=$(subst .c,.o, $(SRCS))

all: sph

%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) -o $@ $^

sph: $(OBJS)
	$(CC) $(CFLAGS) $(LDLIBS) -o $(BIN_DIR)/$@ $^

clean:
	$(RM) $(OBJS)
	$(RM) sph
