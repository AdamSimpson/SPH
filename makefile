CC=mpicc -cc=g++
CLIBS= -lglfw3 -lGLEW -framework OpenGL -framework Cocoa -framework IOkit -framework CoreVideo -L/usr/local/lib -lfreetype
CINCLUDES= -I/usr/local/include/freetype2 -I/Developer/NVIDIA/CUDA-6.5/include/ -I/Users/atj/Downloads -I/Users/atj/SPH
CFLAGS= -DGLFW -O3 -ffast-math -lm -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP -lstdc++
CPPFLAGS = -DGLFW -O3 -ffast-math -lm -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP -DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP

all:
	mkdir -p bin
	cd src && $(CC) $(CINCLUDES) $(CFLAGS) $(CLIBS) ogl_utils.cpp particles_gl.cpp container_gl.cpp camera.cpp mover_gl.cpp font_gl.cpp lodepng.c renderer.cpp gl.cpp tunable_parameters.cpp setup.cpp hash_sort.cpp fluid.cpp communication.cpp -o ../bin/sph.out
clean:
	rm -f ./sph.out
	rm -f ./*.o
