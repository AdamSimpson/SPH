#ifndef fluid_fileio_h
#define fluid_fileio_h

#include "fluid.h"
#include "sph.h"

void writeMPI(fluid_particle_t *particles, int fileNum, param_t *params);
void writeFile(fluid_particle_t *particles, int fileNum, param_t *params);

#endif
