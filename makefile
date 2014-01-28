CC=mpicc

LDFLAGS+=-L$(SDKSTAGE)/opt/vc/lib/ -lGLESv2 -lGLEW -lEGL -lopenmaxil -lbcm_host -lvcos -lvchiq_arm -lpthread -lrt -L../libs/ilclient -L../libs/vgfont
INCLUDES+=-I$(SDKSTAGE)/opt/vc/include/ -I$(SDKSTAGE)/opt/vc/include/interface/vcos/pthreads -I$(SDKSTAGE)/opt/vc/include/interface/vmcs_host/linux -I./ -I../libs/ilclient -I../libs/vgfont
CFLAGS= -DGLES -mfloat-abi=hard -mfpu=vfp -O3 -lm -ffast-math

all:
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) egl_utils.c circles_gl.c renderer.c geometry.c hash.c communication.c fluid.c -o bin/sph.out
debug:
	mkdir -p bin
	mpecc -mpilog -g $(CFLAGS) $(INCLUDES) $(LDFLAGS) egl_utils.c circles_gl.c renderer.c geometry.c hash.c communication.c fluid.c -o bin/sph.out
clean:
	rm -f ./bin/sph.out
	rm -f ./*.o

run: copy
	cd $(HOME) && mpirun --host 192.168.3.100,192.168.3.101,192.168.3.102,192.168.3.103 -n 4 ./sph.out && cd SPH 

copy:
	scp ./bin/sph.out pi0:~/
	scp ./bin/sph.out pi1:~/
	scp ./bin/sph.out pi2:~/
	scp ./bin/sph.out pi3:~/
