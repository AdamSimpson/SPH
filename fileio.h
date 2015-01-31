#ifndef fluid_fileio_h
#define fluid_fileio_h

#include "fluid.h"

void writeMPI(fluid_particle *particles, int fileNum, param *params);
void writeFile(fluid_particle *particles, int fileNum, param *params);

#endif
