#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "fluid.h"
#include "hash.h"
#include "geometry.h"
#include "fileio.h"
#include "communication.h"

int main(int argc, char *argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    createMpiTypes();

    param_t params;
    AABB_t water_volume_global;
    AABB_t boundary_global;
    edge_t edges;
    oob_t out_of_bounds;

    params.rank = rank;
    params.nprocs = nprocs;

    params.g = 9.8;
    params.number_steps = 100;
    params.time_step = 0.01;
    params.number_fluid_particles_global = 1000;
    params.rest_density = 1000.0; // water: kg/m^3

    // Boundary box
    boundary_global.min_x = 0.0;
    boundary_global.max_x = 1.0;
    boundary_global.min_y = 0.0;
    boundary_global.max_y = 1.0;
    boundary_global.min_z = 0.0;
    boundary_global.max_z = 1.0;

    // water volume
    water_volume_global.min_x = 0.1;
    water_volume_global.max_x = 0.5;
    water_volume_global.min_y = 0.1;
    water_volume_global.max_y = 0.5;
    water_volume_global.min_z = 0.1;
    water_volume_global.max_z = 0.5;

    // Mass of each particle
    double volume = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y) * (water_volume_global.max_z - water_volume_global.min_z);

    // Initial spacing between particles
    params.spacing_particle = pow(volume/params.number_fluid_particles_global,1.0/3.0);

    // Smoothing radius, h
    params.smoothing_radius = 2.0*params.spacing_particle;

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
    fluid_particle_t *fluid_particles = malloc(bytes);
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
    initParticles(fluid_particles, neighbors, hash, &water_volume_global, start_x, number_particles_x, &edges, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, params.number_fluid_particles_local, params.smoothing_radius);

    // Initial configuration
    int fileNum=0;
    writeMPI(fluid_particles, fileNum++, &params);

    // Main loop
    // In the current form the particles with be re-hashed and halos re-sent for step 0
    int n;
    double start_time, end_time;

    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();
    for(n=0; n<params.number_steps; n++) {

        printf("Rank %d Entering fluid step %d with %d particles\n",rank, n, params.number_fluid_particles_local);

        apply_gravity(fluid_particles, &params);

        // Advance to predicted position
        predict_positions(fluid_particles, &params);

        // Check boundary partitions
        if (n % 10 == 0)
            checkPartition(fluid_particles, &out_of_bounds, &params);

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(fluid_particles, &out_of_bounds, &params);

        startHaloExchange(fluid_particles, &edges, &params);

        hash_fluid(fluid_particles, neighbors, hash, &params);

        finishHaloExchange(fluid_particles, &edges, &params);

        hash_halo(fluid_particles, neighbors, hash, &params);

        int solve_iterations = 4;
        int si;
        for(si=0; si<solve_iterations; si++)
        {
            compute_densities(fluid_particles, &params);

            calculate_lambda(fluid_particles, neighbors, &params);

            update_dp(fluid_particles, neighbors, &params);
        }

        update_velocities(fluid_particles, &params);

        XSPH_viscosity(fluid_particles, neighbors, &params);

        update_positions(fluid_particles, &params);

        if (n % steps_per_frame == 0)
            writeMPI(fluid_particles, fileNum++, &params);

    }
    end_time = MPI_Wtime();
    printf("Rank %d Elapsed seconds: %f, num particles: %d\n", rank, end_time-start_time, params.number_fluid_particles_local);

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

////////////////////////////////////////////////////////////////////////////
// Smoothing Kernels
///////////////////////////////////////////////////////////////////////////

// (h^2 - r^2)^3 normalized in 3D (poly6)
float W(float r, float h)
{
    if(r > h)
        return 0.0f;

    double C = 315.0/(64.0*M_PI*pow((double)h,9.0));
    float W = C*pow((double)(h*h-r*r), 3.0);
    return W;
}

// Gradient (h-r)^3 normalized in 3D (Spikey)
float del_W(float r, float h)
{
    float C = -45.0/(M_PI * pow((double)h, 6.0));
    float del_W = C*(h-r)*(h-r);
    return del_W;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////

void XSPH_viscosity(fluid_particle_t *fluid_particles, neighbor_t* neighbors, param_t *params)
{
    int i,j;
    unsigned int q_index;
    neighbor_t *n;
    float c = params->c;

    float x_diff, y_diff, z_diff, vx_diff, vy_diff, vz_diff, r_mag, w;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        float partial_sum_x = 0.0f;
        float partial_sum_y = 0.0f;
        float partial_sum_z = 0.0f;
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            vx_diff = fluid_particles[q_index].v_x - fluid_particles[i].v_x;
            vy_diff = fluid_particles[q_index].v_y - fluid_particles[i].v_y;
            vz_diff = fluid_particles[q_index].v_z - fluid_particles[i].v_z;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff * z_diff*z_diff);
            w = W(r_mag, params->smoothing_radius);
            partial_sum_x += vx_diff * w;
            partial_sum_y += vy_diff * w;
            partial_sum_z += vz_diff * w;
        }
        partial_sum_x *= c;
        partial_sum_y *= c;
        partial_sum_z *= c;

        fluid_particles[i].v_x += partial_sum_x;
        fluid_particles[i].v_y += partial_sum_y;
        fluid_particles[i].v_z += partial_sum_z;
    }
}

void compute_densities(fluid_particle_t *fluid_particles, param_t *params)
{
    neighbor_t *neighbors = neighbor_grid->neighbors;

    int i,j;
    unsigned int p_index, q_index;
    neighbor_t *n;
    float h = params->tunable_params.smoothing_radius;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        float x_diff, y_diff, z_diff, r_mag, density;
        density = 0.0f;

        // Own contribution to density
        density += W(0.0f, h);

        // Neighbor contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            if(r_mag <= h)
                density += mass*W(r_mag, h);
        }

        // Update particle density
        fluid_particles[i].density = density;
    }
}

void apply_gravity(fluid_particle_t *fluid_particles, param_t *params)
{
    int i;
    unsigned int p_index;
    float dt = params->time_step;
    float g = -params->g;

    for(i=0; i<(params->number_fluid_particles_local); i++) {
        fluid_particles[i].v_y +=g*dt;
    }
}

void update_dp_positions(fluid_particle_t *fluid_particles, param_t *params)
{
    int i;

    for(i=0; i<(params->number_fluid_particles_local); i++) {
        p_index = fluid_particle_indices[i];
        fluid_particles[i].x_star += fluid_particles[i].dp_x;
        fluid_particles[i].y_star += fluid_particles[i].dp_y;
        fluid_particles[i].z_star += fluid_particles[i].dp_z;

        // Enforce boundary conditions
        boundary_conditions(i, fluid_particles, params);
    }
}

void update_positions(fluid_particle_t *fluid_particles, param_t *params)
{
     int i;

     for(i=0; i<(params->number_fluid_particles_local); i++) {
        p_index = fluid_particle_indices[i];
        fluid_particles[i].x = fluid_particles[i].x_star;
        fluid_particles[i].y = fluid_particles[i].y_star;
        fluid_particles[i].z = fluid_particles[i].z_star;
    }
}

void calculate_lambda(fluid_particle_t *fluid_particles, neighbor_t *neighbor_grid, param_t *params)
{
    neighbor_t *neighbors = neighbor_grid->neighbors;

    int i,j;
    unsigned int q_index;
    neighbor_t *n;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        float Ci = fluid_particles[i].density/params->tunable_params.rest_density - 1.0f;

        float sum_C, x_diff, y_diff, z_diff, r_mag, grad, grad_x, grad_y, grad_z;

        sum_C = 0.0f;
        grad_x = 0.0f;
        grad_y = 0.0f;
        grad_z = 0.0f;
        // Add k = i contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            grad = del_W(r_mag, params->smoothing_radius);
            grad_x += grad*x_diff;
            grad_y += grad*y_diff;
            grad_z += grad*z_diff;
           }
           sum_C += grad_x*grad_x + grad_y*grad_y + grad_z*grad_z;

        // Add k =j contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            grad = del_W(r_mag, params->smoothing_radius);
            grad_x = grad*x_diff ;
            grad_y = grad*y_diff;
            grad_z = grad*z_diff;
            sum_C += (grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
        }

        sum_C *= (1.0f/(params->rest_density*params->rest_density));

        float epsilon = 1.0f;
        fluid_particles[i].lambda = -Ci/(sum_C + epsilon);
    }
}

void update_dp(fluid_particle_t *fluid_particles, neighbor_t *neighbor_grid, param_t *params)
{
    neighbor_t *neighbors = fluid_sim->neighbor_grid->neighbors;

    unsigned int q_index;
    neighbor_t *n;
    float x_diff, y_diff, z_diff, dp, r_mag;

    int i,j;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        float dp_x = 0.0f;
        float dp_y = 0.0f;
        float dp_z = 0.0f;
        float s_corr;
        float k = params->k;
        float dq = params->dq;
        float Wdq = W(dq, params->smoothing_radius);

        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->fluid_neighbors[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            s_corr = -k*(powf(W(r_mag, params->tunable_params.smoothing_radius)/Wdq, 4.0f));
            dp = (fluid_particles[i].lambda + fluid_particles[q_index].lambda + s_corr)*del_W(r_mag, params->smoothing_radius);
            dp_x += dp*x_diff;
            dp_y += dp*y_diff;
            dp_z += dp*z_diff;
        }
        fluid_particles[i].dp_x = dp_x/params->rest_density;
        fluid_particles[i].dp_y = dp_y/params->rest_density;
        fluid_particles[i].dp_z = dp_z/params->rest_density;
    }
}

// Identify out of bounds particles and send them to appropriate rank
void identify_oob_particles(fluid_particle_t *fluid_particles, oob_t *out_of_bounds, param_t *params)
{
    int i;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        // Set OOB particle indices and update number
        if (fluid_particles[i].x < params->tunable_params.node_start_x)
            out_of_bounds->oob_index_indices_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (fluid_particles[i].x > params->tunable_params.node_end_x)
            out_of_bounds->oob_index_indices_right[out_of_bounds->number_oob_particles_right++] = i;
    }

   // Transfer particles that have left the processor bounds
   transfer_OOB_particles(fluid_particles, out_of_bounds);
}

// Predict position
void predict_positions(fluid_particle_t *fluid_particles, param_t *params)
{
    int i;
    float dt = params->time_step;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particles[i].x_star = fluid_particles[i].x + (fluid_particles[i].v_x * dt);
        fluid_particles[i].y_star = fluid_particles[i].y + (fluid_particles[i].v_y * dt);
        fluid_particles[i].z_star = fluid_particles[i].z + (fluid_particles[i].v_z * dt);

        // Enforce boundary conditions
        boundary_conditions(i, fluid_particles);
    }
}

void check_velocity(float *v_x, float *v_y, float *v_z)
{
    const float v_max = 20.0f;

    if(*v_x > v_max)
        *v_x = v_max;
    else if(*v_x < -v_max)
        *v_x = -v_max;
    if(*v_y > v_max)
        *v_y = v_max;
    else if(*v_y < -v_max)
        *v_y = -v_max;
    if(*v_z > v_max)
        *v_z = v_max;
    else if(*v_z < -v_max)
        *v_z = -v_max;

}

// Update particle position and check boundary
void update_velocities(fluid_particle_t *fluid_particles, param_t *params)
{
    float dt = params->time_step;
    float v_x, v_y, v_z;

    // Update local and halo particles, update halo so that XSPH visc. is correct
    for(i=0; i<params->number_fluid_particles_local + params->number_halo_particles; i++) {
        v_x = (fluid_particles[i].x_star - fluid_particles[i].x)/dt;
        v_y = (fluid_particles[i].y_star - fluid_particles[i].y)/dt;
        v_z = (fluid_particles[i].z_star - fluid_particles[i].z)/dt;

        check_velocity(&v_x, &v_y, &v_z);

        fluid_particles[i].v_x = v_x;
        fluid_particles[i].v_y = v_y;
        fluid_particles[i].v_z = v_z;
    }
}

// Assume AABB with min point being axis origin
void boundary_conditions(fluid_particle_t *fluid_particles, unsigned int i, AABB_t *boudnary)
{

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(fluid_particles[i].x_star  < boundary->min_x) {
        fluid_particles[i].x_star = boundary->min_x;
    }
    else if(fluid_particles[i].x_star  > boundary->max_x){
        fluid_particles[i].x_star = boundary->max_x-0.001f;
    }
    if(fluid_particles[i].y_star  <  boundary->min_y) {
        fluid_particles[i].y_star = boundary->min_y;
    }
    else if(fluid_particles[i].y_star  > boundary->max_y){
        fluid_particles[i].y_star = boundary->max_y-0.001f;
    }
    if(fluid_particles[i].z_star  <  boundary->min_z) {
        fluid_particles[i].z_star = boundary->min_z;
    }
    else if(fluid_particles[i].z_star  > boundary->max_z){
        fluid_particles[i].z_star = boundary->max_z-0.001f;
    }
}

// Initialize particles
void initParticles(fluid_particle_t *fluid_particles,
                neighbor *neighbors, n_bucket *hash, AABB* water,
                int start_x, int number_particles_x, edge *edges, param* params)
{
    int i;

    // Create fluid volume
    constructFluidVolume(fluid_particles, water, start_x, number_particles_x, edges, params);

    // Initialize particle values
    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particles[i].v_x = 0.0;
        fluid_particles[i].v_y = 0.0;
        fluid_particles[i].v_z = 0.0;
        fluid_particles[i].density = params->rest_density;
    }

    // Send halo particles
    startHaloExchange(fluid_particles, edges, params);

    //Generate neighbor hash
    hash_fluid(fluid_particles, neighbors, hash, params);
    finishHaloExchange(fluid_particles, edges, params);
    hash_halo(fluid_particles, neighbors, hash, params);
}
