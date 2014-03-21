#ifndef fluid_communication_h
#define fluid_communication_h

typedef struct EDGE_T edge_t;
typedef struct OOB_T oob_t;

#include "fluid.h"
#include "mpi.h"

// MPI globals
MPI_Datatype Particletype;
MPI_Datatype TunableParamtype;
MPI_Datatype LeftEdgetype;
MPI_Datatype RightEdgetype;
MPI_Comm MPI_COMM_COMPUTE;
MPI_Group group_world;
MPI_Group group_compute;
MPI_Group group_render;

// Particles that are within 2*h distance of node edge
struct EDGE_T {
    int max_edge_particles;
    fluid_particle **edge_pointers_left;
    fluid_particle **edge_pointers_right;
    int number_edge_particles_left;
    int number_edge_particles_right;
    MPI_Request reqs[4];
};

// Particles that have left the node
struct OOB_T {
    int max_oob_particles;
    int *oob_pointer_indicies_left; // Indicies in particle pointer array for particles traveling left
    int *oob_pointer_indicies_right;
    int number_oob_particles_left;
    int number_oob_particles_right;
    int *vacant_indicies; // Indicies in particle array that are vacant
    int number_vacancies;
};

void createMpiTypes();
void create_communicators();
void freeMpiTypes();
void startHaloExchange(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,  edge_t *edges, param *params);
void finishHaloExchange(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,  edge_t *edges, param *params);
void transferOOBParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, oob_t *out_of_bounds, param *params);

#endif
