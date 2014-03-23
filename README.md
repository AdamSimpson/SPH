TinySPH
===

TinySPH is a parallel 2D Smoothed Particle Hydrodynamics(SPH) code, designed to run in real time on a Raspberry Pi cluster. This application was designed to demonstrate distributed computing concepts and while it works best in a distributed environment it is possible to run on a shared memory system such as a single multicore processor. The code is in development but should compile and run on Raspbian, various flavors of linux(untested), and OS X. The code makes use of MPI for distributed computing and requires at least two mpi processes(1 render, 1 compute).

Install
===

Several prerequisites are required before compiling the code. In a Linux environment, such as Raspian, these may be obtained using your distros package management system. On Macintosh it is recomended that Homebrew be used: http://brew.sh

Macintosh
---
Once Homebrew has been installed The following commands may be run in Terminal.app to setup the prerequisties.

$ brew install mpich
$ brew install glfw
$ brew install glew
$ brew install freetype

Raspberry Pi
---
On the RaspberryPi the following packages must be installed, all of which are availble through apt-get.

$ sudo apt-get install mpich
$ sudo apt-get install glew

Compile and run
===
Once the prerequisites have been setup TitanSPH can be compiled.

Macintosh
---
$ make -f makefile_macos

To run on a 4 core single socket computer:

$ mpirun -n 4 ./bin/SPH

Raspbery Pi
---
$ make

To run on TinyTitan

$ make run

