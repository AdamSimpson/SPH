#include "mpi.h"
#include "renderer.hpp"
#include "fluid.h"


int main(int argc, char *argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank;

    // Rank in world space
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    create_communicators();

    create_MPI_types();

    // Rank 0 is the render node, otherwise a simulation node
    if(rank == 0) {
        Renderer renderer;
        renderer.start_renderer();
    }
    else
        start_simulation();

    free_MPI_types();
    MPI_Finalize();
    return 0;
}
