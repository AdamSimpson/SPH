#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

#include "hash.h"
#include "renderer.h"
#include "geometry.h"
#include "fluid.h"
#include "communication.h"

int main(int argc, char *argv[])
{
     // Initialize MPI
    MPI_Init(&argc, &argv);
    int rank;

    // Rank in world space
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    create_communicators();

    createMpiTypes();

    // Rank 0 is the render node, otherwise a simulation node
    if(rank == 0)
	start_renderer();
    else
	start_simulation();
}

void start_simulation()
{
    int rank, nprocs;

    MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
    MPI_Comm_size(MPI_COMM_COMPUTE, &nprocs);

    printf("compute rank: %d, num compute procs: %d \n",rank, nprocs);

    param params;
    AABB water_volume_global;
    AABB boundary_global;
    edge edges;
    oob out_of_bounds;

    params.rank = rank;
    params.nprocs = nprocs;

    params.g = 1.0;
    params.time_step = 0.03;
    // The number of particles used may differ slightly
    params.number_fluid_particles_global = 5000;
    params.rest_density = 10.0;

    // Boundary box
    boundary_global.min_x = 0.0;
    boundary_global.max_x = 20.0;
    boundary_global.min_y = 0.0;
    boundary_global.max_y = 10.0;

    // water volume
    water_volume_global.min_x = 0.0;
    water_volume_global.max_x = 10.0;
    water_volume_global.min_y = 0.0;
    water_volume_global.max_y = 10.0;

    // Fluid area in initial configuration
    double area = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y);

    // Initial spacing between particles
    params.spacing_particle = pow(area/params.number_fluid_particles_global,1.0/2.0);

    // Smoothing radius, h
    params.smoothing_radius = 2.0*params.spacing_particle;

    printf("smoothing radius: %f\n", params.smoothing_radius);

    int start_x;  // where in x direction this nodes particles start
    int number_particles_x; // number of particles in x direction for this node

    // Divide problem set amongst nodes
    partitionProblem(&boundary_global, &water_volume_global, &start_x, &number_particles_x, &params);

    // Set local/global number of particles to allocate
    setParticleNumbers(&boundary_global, &water_volume_global, &edges, &out_of_bounds, number_particles_x, &params);

    long long total_bytes = 0;
    size_t bytes;
    // Allocate fluid particles array
    bytes = params.max_fluid_particles_local * sizeof(fluid_particle);
    total_bytes+=bytes;
    fluid_particle *fluid_particles = malloc(bytes);
    if(fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");

    // Allocate (x,y) coordinate array
    bytes = 2 * params.max_fluid_particles_local * sizeof(float);
    total_bytes+=bytes;
    float *fluid_particle_coords = malloc(bytes);
    if(fluid_particle_coords == NULL)
        printf("Could not allocate fluid_particle coords\n");

    // Allocate pointer array used to traverse non vacant particles
    bytes = params.max_fluid_particles_local * sizeof(fluid_particle*);
    total_bytes+=bytes;
    fluid_particle **fluid_particle_pointers = malloc(bytes);
    if(fluid_particle_pointers == NULL)
        printf("Could not allocate fluid_particle_pointers\n");

    // Allocate neighbors array
    neighbor *neighbors = calloc(params.max_fluid_particles_local, sizeof(neighbor));
    total_bytes+=params.max_fluid_particles_local*sizeof(neighbor);
    if(neighbors == NULL)
        printf("Could not allocate neighbors\n");

    // UNIFORM GRID HASH
    params.grid_size_x = ceil((boundary_global.max_x - boundary_global.min_x) / params.smoothing_radius);
    params.grid_size_y = ceil((boundary_global.max_y - boundary_global.min_y) / params.smoothing_radius);
    unsigned int grid_size = params.grid_size_x * params.grid_size_y;
    params.length_hash = grid_size;
    n_bucket* hash = calloc(params.length_hash, sizeof(n_bucket));
    if(hash == NULL)
        printf("Could not allocate hash\n");

    // Allocate edge index arrays
    edges.edge_pointers_left = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    edges.edge_pointers_right = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    // Allocate out of bound index arrays
    out_of_bounds.oob_pointer_indicies_left = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.oob_pointer_indicies_right = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.vacant_indicies = malloc(2*out_of_bounds.max_oob_particles * sizeof(int));

    printf("bytes allocated: %lld\n", total_bytes);

    // Initialize particles
    initParticles(fluid_particle_pointers, fluid_particles, neighbors, hash, &water_volume_global, start_x, number_particles_x, &edges, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, params.number_fluid_particles_local, params.smoothing_radius);

    // Send intiial paramaters to render node
    param *null_param = NULL;
    int *null_recvcnts = NULL;
    int *null_displs = NULL;
    MPI_Gatherv(&params, 1, Paramtype, null_param, null_recvcnts, null_displs, Paramtype, 0, MPI_COMM_WORLD);

    double start_time, end_time, partition_time;
    fluid_particle *p;
    unsigned int i;
    unsigned int n = 0;
    fluid_particle *null_particle;

    // Main simulation loop
    while(1) {
        // Initialize velocities
	apply_gravity(fluid_particle_pointers, &params);

        // Viscosity impluse
        // This is missing halo particle contribution
        viscosity_impluses(fluid_particle_pointers, neighbors, &params);

        // Advance to predicted position and set OOB particles
        predict_positions(fluid_particle_pointers, &out_of_bounds, &boundary_global, &params);

        // Transfer particles that have left the processor bounds
        transferOOBParticles(fluid_particle_pointers, fluid_particles, &out_of_bounds, &params);

        // Send compute parameters to render node
        MPI_Gatherv(&params, 1, Paramtype, null_param, null_recvcnts, null_displs, Paramtype, 0, MPI_COMM_WORLD);
        // Receive updated paramaters from render nodes
        MPI_Scatterv(null_param, 0, null_displs, Paramtype, &params, 1, Paramtype, 0,  MPI_COMM_WORLD);

        // Hash the non halo regions
	// This will update the densities so when the halo is exchanged the halo particles are up to date
	// This works well on the raspi's but destroys communication/computation overlap
        hash_fluid(fluid_particle_pointers, neighbors, hash, &params);

	// Exchange halo particles
        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);
        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

	// Every 10 steps start the timer used to determine if the partition needs to be modified
	if(n % 10 == 0)
	    partition_time = MPI_Wtime();

	// Add the halo particles to neighbor buckets
	// Also update density
        hash_halo(fluid_particle_pointers, neighbors, hash, &params);

        // double density relaxation
	// halo particles will be missing origin contributions to density/pressure
        double_density_relaxation(fluid_particle_pointers, neighbors, &params);

        // update velocity
        updateVelocities(fluid_particle_pointers, &edges, &boundary_global, &params);

	// Check for a balanced particle load between MPI tasks
        if (n % 10 == 0) {
            checkPartition(fluid_particle_pointers, &out_of_bounds, &partition_time, &params);
            // reset loop count
	    n = 0;
        }

        // Pack fluid_particle_coords
	for(i=0; i<params.number_fluid_particles_local; i++) {
	    p = fluid_particle_pointers[i];
	    fluid_particle_coords[i*2] = p->x;
	    fluid_particle_coords[(i*2)+1] = p->y;
        }

        // Send particle positions to be rendered
        MPI_Gatherv(fluid_particle_coords, 2*params.number_fluid_particles_local, MPI_FLOAT,
		    null_particle, null_recvcnts, null_displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
 
       // iterate sim loop counter
	n++;
    }

    // Release memory
    free(fluid_particle_pointers);
    free(fluid_particles);
    free(neighbors);
    free(hash);

    // Close MPI
    freeMpiTypes();
    MPI_Finalize();

}

// This should go into the hash, perhaps with the viscocity?
void apply_gravity(fluid_particle **fluid_particle_pointers, param *params)
{
    int i;
    fluid_particle *p;
    double dt = params->time_step;
    double g = -params->g;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        p->v_y += g*dt;

        // Zero out density as well
        p->density = 0.0;
        p->density_near = 0.0;
     }
}

void viscosity_impulse(fluid_particle *p, fluid_particle *q, param* params)
{

    static const double sigma = 0.5;
    static const double beta =  0.1;
    double h = params->smoothing_radius;
    double dt = params->time_step;

    double r, ratio, u, imp, imp_x, imp_y;

    r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y));
    ratio = r/h;

    //Inward radial velocity
    u = ((p->v_x-q->v_x)*(q->x-p->x) + (p->v_y-q->v_y)*(q->y-p->y))/r;
    if(u>0.0)
    {
	imp = dt * (1-ratio)*(sigma * u + beta * u*u);
	imp_x = imp*(q->x-p->x)/r;
	imp_y = imp*(q->y-p->y)/r;
	p->v_x -= imp_x/2.0;
	p->v_y -= imp_y/2.0;
	q->v_x += imp_x/2.0;
	q->v_y += imp_y/2.0;
    }

}

// Add viscosity impluses
void viscosity_impluses(fluid_particle **fluid_particle_pointers, neighbor* neighbors, param *params)
{
    int i, j;
    fluid_particle *p, *q;
    neighbor* n;
    double r,ratio,h,u,dt;
    double imp, imp_x, imp_y;
    static const double sigma = 0.5;
    static const double beta =  0.1;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];
        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
	    viscosity_impulse(p,q,params);
        }
    }
}

// Predict position
void predict_positions(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, AABB *boundary_global, param *params)
{
    int i;
    fluid_particle *p;
    double dt = params->time_step;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
	p->x_prev = p->x;
        p->y_prev = p->y;
	p->x += (p->v_x * dt);
        p->y += (p->v_y * dt);

	// Enforce boundary conditions
        boundaryConditions(p, boundary_global, params);

        // Set OOB particle indicies and update number
        if (p->x < params->node_start_x)
            out_of_bounds->oob_pointer_indicies_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (p->x >= params->node_end_x)
            out_of_bounds->oob_pointer_indicies_right[out_of_bounds->number_oob_particles_right++] = i;
    }
}

void calculate_density(fluid_particle *p, fluid_particle *q, param *params)
{
    double ratio, r;
    r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y));
    ratio = r/params->smoothing_radius;
    if(ratio < 1.0) {
	p->density += (1-ratio)*(1-ratio);
	p->density_near += (1-ratio)*(1-ratio)*(1-ratio);

	q->density += (1-ratio)*(1-ratio);
	q->density_near += (1-ratio)*(1-ratio)*(1-ratio);
    }

}

void double_density_relaxation(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    static const double k = 0.2;
    static const double k_near = 2.0;

    int i, j;
    fluid_particle *p, *q;
    neighbor* n;
    double r,ratio,dt,h,u,dx_x,dx_y,D,D_x,D_y;
    h = params->smoothing_radius;
    dt = params->time_step;

    // Calculate the pressure of all particles, including halo
    for(i=0; i<params->number_fluid_particles_local + params->number_halo_particles; i++) {
        p = fluid_particle_pointers[i];
        // Compute pressure and near pressure
        p->pressure = k * (p->density-params->rest_density);
        p->pressure_near = k_near * p->density_near;
    }

    // Relax particle positions
//    for(i=0; i<params->number_fluid_particles_local; i++) {

    // Iterating through the array in reverse seems to have a noticable effect on stability
    for(i=params->number_fluid_particles_local; i-- > 0; ) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        for(j=0; j<n->number_fluid_neighbors; j++) {

            q = n->fluid_neighbors[j];
            r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y));
	    ratio = r/h;
	    if(ratio < 1.0 && ratio > 0.0) {
		// Updating both neighbor pairs at the same time, slightly different than the paper but quicker
	        // Also the running sum of D for particle p seems to produce more bias/instability so is removed
                D = dt*dt*((p->pressure+q->pressure)*(1.0-ratio) + (p->pressure_near+q->pressure_near)*(1.0-ratio)*(1.0-ratio));
		D_x = D*(q->x-p->x)/r;
                D_y = D*(q->y-p->y)/r;

		// Do not move the halo particles full D
		// Halo particles are missing D from their origin so I believe this is appropriate
		if(q->id < params->number_fluid_particles_local) {
	  	  q->x += D_x;
	          q->y += D_y;
		}	
		else { // Not even sure what i'm doing anymore
                  q->x += D_x/2.0;
                  q->y += D_y/2.0;
		}
		    

		p->x -= D_x;
                p->y -= D_y;
            }
        }

    }
}

void updateVelocity(fluid_particle *p, param *params)
{
    double dt = params->time_step;

    p->v_x = (p->x-p->x_prev)/dt;
    p->v_y = (p->y-p->y_prev)/dt;
}

// Update particle position and check boundary
void updateVelocities(fluid_particle **fluid_particle_pointers, edge *edges, AABB *boundary_global, param *params)
{
    int i;
    fluid_particle *p;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        boundaryConditions(p, boundary_global, params);
        updateVelocity(p, params);

    }
}

void collisionImpulse(fluid_particle *p, int norm_x, int norm_y, param *params)
{
    // Boundary friction
    const static double mu = 0.0;
    double dt = params->time_step;

    double v_dot_n, vx_norm, vy_norm, vx_tan, vy_tan, I_x, I_y;

    // v.n
    v_dot_n = p->v_x*norm_x + p->v_y*norm_y;

    // Velocity normal to surface
    vx_norm = v_dot_n*norm_x;
    vy_norm = v_dot_n*norm_y;

    // Velocity tangential to surface
    vx_tan = p->v_x - vx_norm;
    vy_tan = p->v_y - vy_norm;

    // Impulse
    I_x = vx_norm - mu*vx_tan;
    I_y = vy_norm - mu*vy_tan;

    p->x += I_x * dt;
    p->y += I_y * dt;
}

// Assume AABB with min point being axis origin
void boundaryConditions(fluid_particle *p, AABB *boundary, param *params)
{
    // Update velocity
    if(p->x < boundary->min_x) {
	collisionImpulse(p,1,0,params);
    }
    else if(p->x > boundary->max_x){
        collisionImpulse(p, -1, 0, params);
    }
    if(p->y < boundary->min_y) {
        collisionImpulse(p,0,1,params);
    }
    else if(p->y > boundary->max_y){
        collisionImpulse(p,0,-1,params);
    }

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash won't pick it up
    // as the particle will in the 'next' after last x direction bin
    if(p->x < boundary->min_x) {
//        p->x_prev = p->x;
        p->x = boundary->min_x;
    }
    else if(p->x > boundary->max_x){
//        p->x_prev = p->x;
        p->x = boundary->max_x-0.001;
    }
    if(p->y < boundary->min_y) {
//        p->y_prev = p->y;
        p->y = boundary->min_y;
    }
    else if(p->y > boundary->max_y){
//        p->y_prev = p->y;
        p->y = boundary->max_y-0.001;
    }

}

// Initialize particles
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                neighbor *neighbors, n_bucket *hash, AABB* water, int start_x, int number_particles_x, edge *edges, param* params)
{
    int i;
    fluid_particle *p;

    // Create fluid volume
    constructFluidVolume(fluid_particle_pointers, fluid_particles, water, start_x, number_particles_x, edges, params);

    // NULL out unused fluid pointers
    for(i=params->number_fluid_particles_local; i<params->max_fluid_particles_local; i++)
        fluid_particle_pointers[i] = NULL;

    // Initialize particle values
    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particle_pointers[i]->a_x = 0.0;
        fluid_particle_pointers[i]->a_y = 0.0;
        fluid_particle_pointers[i]->v_x = 0.0;
        fluid_particle_pointers[i]->v_y = 0.0;
    }

}
