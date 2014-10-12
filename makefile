CC=mpicc
CLIBS= -lglfw3 -lGLEW -framework OpenGL -framework Cocoa -framework IOkit -framework CoreVideo -L/usr/local/lib -lfreetype
CINCLUDES= -I/usr/local/include/freetype2
CFLAGS= -DGLFW -O3 -ffast-math -lm

all:
	mkdir -p bin
	cd src && $(CC) $(CINCLUDES) $(CFLAGS) $(CLIBS) ogl_utils.c dividers_gl.c particles_gl.c liquid_gl.c mover_gl.c font_gl.c lodepng.c renderer.c glfw_utils.c image_gl.c cursor_gl.c background_gl.c controls.c setup.c hash.c communication.c fluid.c -o ../bin/sph.out
clean:
	rm -f ./sph.out
	rm -f ./*.o
