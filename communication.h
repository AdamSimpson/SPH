#ifndef fluid_communication_h
#define fluid_communication_h

typedef struct EDGES edge_t;
typedef struct OOB oob_t;

#include "fluid.h"
#include "mpi.h"

MPI_Datatype Particletype;

// Particles that are within halo width of node edge
struct EDGES {
    int max_edge_particles;
    int *edge_indices_left;
    int *edge_indices_right;
    int number_edge_particles_left;
    int number_edge_particles_right;
    MPI_Request reqs[4];
};

// Particles that have left the node
struct OOB {
    int max_oob_particles;
    int *oob_indices_left; // Indicies in particle pointer array for particles traveling left
    int *oob_indices_right;
    int number_oob_particles_left;
    int number_oob_particles_right;
};

void createMpiTypes();
void freeMpiTypes();
void transferHalos(fluid_particle_t *fluid_particles, edge_t *edges, param_t *params);
void transferOOBParticles(fluid_particle_t *fluid_particles, oob_t *out_of_bounds, param_t *params);
void startHaloExchange(fluid_particle_t *fluid_particles, edge_t *edges, param_t *params);
void finishHaloExchange(fluid_particle_t *fluid_particles, edge_t *edges, param_t *params);
void update_halo_lambdas(fluid_particle_t *fluid_particles,  edge_t *edges, param_t *params);
void update_halo_positions(fluid_particle_t *fluid_particles,  edge_t *edges, param_t *params);

#endif
