# TinySPH
TinySPH is a parallel 2D Smoothed Particle Hydrodynamics(SPH) code, designed to run in real time on the Oak Ridge Leadership Computing Facilities "Tiny Titan" Raspberry Pi cluster. This application was designed to demonstrate distributed computing concepts and while it works best in a distributed environment it is possible to run on a shared memory system such as a single multicore processor. Although the code is designed for Tiny Titan it should, perhaps with some small build tweaks, compile and run on various flavors of Linux and OS X. The code makes use of MPI for distributed computing and requires at least two mpi processes(1 render, 1 compute).

If you find this code useful, find a bug, or use this code in an interesting way i'd love to hear about it, drop me a line at simpsonab@ornl.gov and let me know!

The screenshot below shows TinySPH running on a four core MacBook Pro, running with 3 compute MPI processes and 1 MPI render process. The color of the particle indicates which processor core is responsable for it's computation. TinySPH includes a simple number of particle based load balancing scheme.

![alt text](https://raw.githubusercontent.com/AdamSimpson/SPH/master/SPH_Screenshot.png "SPH Screenshot")

## Install

Several prerequisites are required before compiling the code. In a Linux environment, such as Raspian, these may be obtained using your distros package management system. On Macintosh it is recomended that [Homebrew](http://brew.sh) be used.

### Macintosh

It is assumed that the XCode toolchain has been installed, this is freely available from [Apple](https://developer.apple.com/xcode/downloads/) . Once Homebrew has been installed The following commands may be run in Terminal.app to setup the prerequisties.

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
The input controls are set in `GLFW_utils.c` and `EGL_utils.c` for GLFW and Raspberry Pi platforms respectively. The Pi's controls are based upon using an XBox controller to handle input.

Although subject to change generally keys should operate as follows:

* The escape key should exit the simulation, if all else failse control-c will do it as well

* Arrow keys should change the paramaters in the top left of the screen

* The mouse controlls the mover sphere

* A,B,X,Y are fluid parameter presets

* on Mac [ ] control the number of processes while on the Pi it is page up and page down.

If the keyboard input for the RaspberyPi doesn't work you may need to correctly set `/dev/input/event#` in `get_key_press()` in `egl_util.c` 

## Code
If you wish the modify the code here are a few things to keep in mind

* `main()` lives in `fluid.c` ~line 37
* The main computation loop lives in `fluid.c` ~lines 236
* Initial parameters are largely set in `fluid.c` starting ~line 78
* The shaders directory contains OpenGL and OpenGL ES 2.0 shaders
* files with suffix \_gl control OpenGL rendering

## Algorithm
The SPH algorithm is based upon the work of [Clavet et al.](http://www.ligum.umontreal.ca/Clavet-2005-PVFS/pvfs.pdf). To run in real time on the Raspberry Pi a large timestep was neccessary as communication is extremely expensive. Several modifcations have been made to make the algorithm work on the RaspberryPi.

* 2D was used instead of 3D.
* The algorithm was parallelized, the prediction-relaxation algorithm does not have a trivial parallelization scheme.
* The algorithm was completely symmetrized to reduce computation/hash cost.
* Plasticity is not used.
