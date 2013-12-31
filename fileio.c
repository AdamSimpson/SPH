////////////////////////////////////////////////
// File input/output functions
////////////////////////////////////////////////

#include <stdio.h>
#include "fileio.h"
#include "fluid.h"

// Write fluid particle data to file
void writeFile(fluid_particle *particles, int fileNum, param *params)
{
    fluid_particle *p;
    FILE *fp ;
    int i;
    char name[64];
    sprintf(name, "/tmp/work/atj/sim-%d.csv", fileNum);
    fp = fopen ( name,"w" );
    if (!fp) {
        printf("ERROR: error opening file\n");
        exit(1);
    }
    for(i=0; i<params->number_fluid_particles; i++) {
        p = &particles[i];
        fprintf(fp,"%f,%f,%f\n",p->x,p->y,p->z);
    }
    fclose(fp);
    printf("wrote file: %s\n", name);
}

