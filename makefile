CC=mpicc

LDFLAGS+=-L$(SDKSTAGE)/opt/vc/lib/ -lGLESv2 -lGLEW -lEGL -lopenmaxil -lbcm_host -lvcos -lvchiq_arm -lpthread -lrt -L../libs/ilclient -L../libs/vgfont -lfreetype
INCLUDES+=-I$(SDKSTAGE)/opt/vc/include/ -I$(SDKSTAGE)/opt/vc/include/interface/vcos/pthreads -I$(SDKSTAGE)/opt/vc/include/interface/vmcs_host/linux -I./ -I../libs/ilclient -I../libs/vgfont -I/usr/include/freetype2
CFLAGS= -DRASPI -mfloat-abi=hard -mfpu=vfp -O3 -lm -ffast-math -g

all:
	mkdir -p bin
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) ogl_utils.c egl_utils.c rgb_light.c dividers_gl.c lodepng.c background_gl.c font_gl.c particles_gl.c mover_gl.c controls.c renderer.c geometry.c hash.c communication.c fluid.c -o bin/sph.out

debug:
	mkdir -p bin
	mpecc -mpilog -Wall -g -DRASPI -DLIGHT $(INCLUDES) $(LDFLAGS) ogl_utils.c egl_utils.c rgb_light.c dividers_gl.c lodepng.c background_gl.c font_gl.c particles_gl.c mover_gl.c controls.c renderer.c geometry.c hash.c communication.c fluid.c -o bin/sph.out


clean:
	rm -f ./bin/sph.out
	rm -f ./*.o

run: copy
	cd $(HOME) ; mpirun -f ~/pi_mpihostsfile -n 12 ~/sph.out ; cd $(HOME)/SPH

copy:
	scp ./bin/sph.out pi1:~/
	scp ./bin/sph.out pi2:~/
	scp ./bin/sph.out pi3:~/
	scp ./bin/sph.out pi4:~/
	scp ./bin/sph.out pi5:~/
	scp ./bin/sph.out pi6:~/
	scp ./bin/sph.out pi7:~/
	scp ./bin/sph.out pi8:~/
	scp ./bin/sph.out pi9:~/
	scp ./bin/sph.out pi10:~/
	scp ./bin/sph.out pi11:~/
	scp ./bin/sph.out pi12:~/
