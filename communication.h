#ifndef fluid_communication_h
#define fluid_communication_h

typedef struct EDGES edge;
typedef struct OOB oob;

#include "fluid.h"
#include "mpi.h"

MPI_Datatype Particletype;
MPI_Datatype LeftEdgetype;
MPI_Datatype RightEdgetype;

// Particles that are within 2*h distance of node edge
struct EDGES {
    int max_edge_particles;
    fluid_particle **edge_pointers_left;
    fluid_particle **edge_pointers_right;
    int number_edge_particles_left;
    int number_edge_particles_right;
    MPI_Request reqs[4];
};

// Particles that have left the node
struct OOB {
    int max_oob_particles;
    int *oob_pointer_indicies_left; // Indicies in particle pointer array for particles traveling left
    int *oob_pointer_indicies_right;
    int number_oob_particles_left;
    int number_oob_particles_right;
    int *vacant_indicies; // Indicies in particle array that are vacant
    int number_vacancies;
};

void createMpiTypes();
void freeMpiTypes();
void transferHalos(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, edge *edges, param *params);
void transferOOBParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, oob *out_of_bounds, param *params);
void startHaloExchange(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,  edge *edges, param *params);
void finishHaloExchange(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,  edge *edges, param *params);

#endif
