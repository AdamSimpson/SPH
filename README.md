# TinySPH

TinySPH is a parallel 2D Smoothed Particle Hydrodynamics(SPH) code, designed to run in real time on the Oak Ridge Leadership Computing Facilities "Tiny Titan" Raspberry Pi cluster. This application was designed to demonstrate distributed computing concepts and while it works best in a distributed environment it is possible to run on a shared memory system such as a single multicore processor. Although the code is designed for Tiny Titan it should, perhaps with some small build tweaks, compile and run on various flavors of linux and OS X. The code makes use of MPI for distributed computing and requires at least two mpi processes(1 render, 1 compute).

## Install

Several prerequisites are required before compiling the code. In a Linux environment, such as Raspian, these may be obtained using your distros package management system. On Macintosh it is recomended that Homebrew be used: http://brew.sh

### Macintosh

It is assumed that the XCode toolchain has been installed, this is freely available from (https://developer.apple.com/xcode/downloads/ "Apple") . Once Homebrew has been installed The following commands may be run in Terminal.app to setup the prerequisties.

    $ brew install mpich
    $ brew install glfw
    $ brew install glew
    $ brew install freetype

### Raspberry Pi

On the RaspberryPi the following packages must be installed, all of which are availble through apt-get.

    $ sudo apt-get install mpich
    $ sudo apt-get install glew

## Compile and run
Once the prerequisites have been installed TitanSPH can be compiled and run.

### Macintosh

To compile

    $ make -f makefile_macos

To run on a 4 core single socket computer:

    $ mpirun -n 4 ./bin/SPH

### Raspbery Pi
To Compile

    $ make

To run on TinyTitan

    $ make run

## Controls
The controls are listed in GLFW_utils.c and EGL_utils.c for GLFW and Raspberry Pi platforms respectively. The Pi's controls are based upon using an XBox controller to handle input.

Although subject to change generally keys should operate as follows:

* Arrow keys should change the paramaters in the top left of the screen

* The mouse controlls the mover sphere

* A,B,X,Y are fluid parameter presets

* on Mac [ ] control the number of processes while on the Pi it is page up and page down.
