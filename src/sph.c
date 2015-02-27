#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "communication.h"
#include "fluid.h"
#include "geometry.h"
#include "fileio.h"
#include "sph.h"

int main(int argc, char *argv[])
{
    init_communication(argc, argv);

    param_t params;
    AABB_t water_volume_global;
    AABB_t boundary_global;
    edge_t edges;
    oob_t out_of_bounds;

    params.rank = get_rank();
    params.nprocs = get_num_procs(0);
    params.g = 9.8;
    params.number_steps = 2000;
    params.time_step = 1.0/60.0;
    params.c = 0.01;
    params.k = 0.1;
    params.number_fluid_particles_global = 65536*2;

    boundary_global.min_x = 0.0;
    boundary_global.max_x = 100.0;
    boundary_global.min_y = 0.0;
    boundary_global.max_y = 80.0;
    boundary_global.min_z = 0.0;
    boundary_global.max_z = 30.0;

    water_volume_global.min_x = 0.1;
    water_volume_global.max_x = boundary_global.max_x - 20.0;
    water_volume_global.min_y = 0.1;
    water_volume_global.max_y = boundary_global.max_y - 30.0;
    water_volume_global.min_z = 0.1;
    water_volume_global.max_z = boundary_global.max_z - 10.0;

    // Cubed volume
    double volume = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y) * (water_volume_global.max_z - water_volume_global.min_z);

    // Initial spacing between particles
    float spacing_particle = pow(volume/params.number_fluid_particles_global,1.0/3.0);

    // Let mass of each particle equal 1
    params.rest_density = params.number_fluid_particles_global/volume;
    printf("rest density: %f\n", params.rest_density);

    // Smoothing radius, h
    params.smoothing_radius = 2.0*spacing_particle;
    params.dq = 0.1*params.smoothing_radius;

    printf("smoothing radius: %f\n", params.smoothing_radius);

    // Number of steps before frame needs to be written for 30 fps
    int steps_per_frame = (int)(1.0/(params.time_step*30.0));

    int start_x;  // where in x direction this nodes particles start
    int number_particles_x; // number of particles in x direction for this node
    // Divide problem set amongst nodes
    partitionProblem(&boundary_global, &water_volume_global, &start_x, &number_particles_x, &params);
    // Set local/global number of particles to allocate
    setParticleNumbers(&boundary_global, &water_volume_global, &edges, &out_of_bounds, number_particles_x, &params);

    long long total_bytes = 0;
    size_t bytes;
    // Allocate fluid particles array
    bytes = params.max_fluid_particles_local * sizeof(fluid_particle_t);
    total_bytes+=bytes;
    fluid_particle_t *fluid_particles = calloc(params.max_fluid_particles_local, sizeof(fluid_particle_t));
    if(fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");

    // Allocate neighbors array
    neighbor_t *neighbors = calloc(params.max_fluid_particles_local, sizeof(neighbor_t));
    total_bytes+=params.max_fluid_particles_local*sizeof(neighbor_t);
    if(neighbors == NULL)
        printf("Could not allocate neighbors\n");
    /////////////////////
    // UNIFORM GRID HASH
    /////////////////////
    // +1 added because range begins at 0
    params.grid_size_x = ceil((boundary_global.max_x - boundary_global.min_x) / params.smoothing_radius) + 1;
    params.grid_size_y = ceil((boundary_global.max_y - boundary_global.min_y) / params.smoothing_radius) + 1;
    params.grid_size_z = ceil((boundary_global.max_z - boundary_global.min_z) / params.smoothing_radius) + 1;
    unsigned int grid_size = params.grid_size_x * params.grid_size_y * params.grid_size_z;
    params.length_hash = grid_size;
    bucket_t* hash = calloc(params.length_hash, sizeof(bucket_t));
    if(hash == NULL)
        printf("Could not allocate hash\n");

    // Allocate edge index arrays
    edges.edge_indices_left = malloc(edges.max_edge_particles * sizeof(int));
    edges.edge_indices_right = malloc(edges.max_edge_particles * sizeof(int));
    // Allocate out of bound index arrays
    out_of_bounds.oob_indices_left = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.oob_indices_right = malloc(out_of_bounds.max_oob_particles * sizeof(int));

    printf("gigabytes allocated: %lld\n", total_bytes/1073741824);

    // Initialize particles
    initParticles(fluid_particles, neighbors, hash, &water_volume_global, &boundary_global, start_x, number_particles_x, &edges, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", params.rank, params.number_fluid_particles_local, params.smoothing_radius);

    // Initial configuration
    int fileNum=0;
    writeMPI(fluid_particles, fileNum++, &params);

    // Main loop
    int n;
    double start_time, end_time;

    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();
    for(n=0; n<params.number_steps; n++) {

        printf("Rank %d Entering fluid step %d with %d particles\n",params.rank, n, params.number_fluid_particles_local);

        apply_gravity(fluid_particles, &params);

        // Advance to predicted position
        predict_positions(fluid_particles, &boundary_global, &params);

        // Check boundary partitions
        if (n % 10 == 0)
            checkPartition(fluid_particles, &out_of_bounds, &params);

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(fluid_particles, &out_of_bounds, &params);

        startHaloExchange(fluid_particles, &edges, &params);

        hash_fluid(fluid_particles, neighbors, hash, &boundary_global, &params);

        finishHaloExchange(fluid_particles, &edges, &params);

        hash_halo(fluid_particles, neighbors, hash, &boundary_global, &params);

        int solve_iterations = 4;
        int si;
        for(si=0; si<solve_iterations; si++)
        {
            compute_densities(fluid_particles, neighbors, &params);

            calculate_lambda(fluid_particles, neighbors, &params);

            update_halo_lambdas(fluid_particles, &edges, &params);

            update_dp(fluid_particles, neighbors, &params);

            update_dp_positions(fluid_particles, &boundary_global, &params);

            update_halo_positions(fluid_particles, &edges, &params);
        }

        update_velocities(fluid_particles, &params);

        XSPH_viscosity(fluid_particles, neighbors, &params);

        vorticity_confinement(fluid_particles, neighbors, &params);

        update_positions(fluid_particles, &params);

        if (n % steps_per_frame == 0)
            writeMPI(fluid_particles, fileNum++, &params);

    }
    end_time = MPI_Wtime();
    printf("Rank %d Elapsed seconds: %f, num particles: %d\n", params.rank, end_time-start_time, params.number_fluid_particles_local);

    // Release memory
    free(fluid_particles);
    free(neighbors);
    free(hash);
    free(edges.edge_indices_left);
    free(edges.edge_indices_right);
    free(out_of_bounds.oob_indices_left);
    free(out_of_bounds.oob_indices_right);

    // Close MPI
    freeMpiTypes();
    MPI_Finalize();

    return 0;
}
