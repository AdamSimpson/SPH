#include "mpi.h"
#include "renderer.hpp"
#include "fluid.h"


int main(int argc, char *argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank, size;

    // Rank in world space
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Total size in world space
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    create_communicators();

    create_MPI_types();

    // Rank 0 is the render node, otherwise a simulation node
    if(rank == 0) {
        int num_compute_procs = size-1;
        Renderer renderer(num_compute_procs);
        renderer.start_renderer();
    }
    else
        start_simulation();

    free_MPI_types();
    MPI_Finalize();
    return 0;
}
