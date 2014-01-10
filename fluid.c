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

    param params;
    AABB water_volume_global;
    AABB boundary_global;
    edge edges;
    oob out_of_bounds;

    params.rank = rank;
    params.nprocs = nprocs;

    params.g = 9.8;
    params.number_steps = 200;
    params.time_step = 0.01;
    params.surface_tension =  0.01;
    params.number_fluid_particles_global = 1000;
    params.rest_density = 100.0; // water: kg/m^3 (This doesn't really make sense for 2D obviously)

    // Boundary box
    boundary_global.min_x = 0.0;
    boundary_global.max_x = 1.0;
    boundary_global.min_y = 0.0;
    boundary_global.max_y = 1.0;

    // water volume
    water_volume_global.min_x = 0.1;
    water_volume_global.max_x = 1.0;
    water_volume_global.min_y = 0.0;
    water_volume_global.max_y = 1.0;

    // Mass of each particle
    // This doesn't really make sense in 2D
    // Try to come up with the mass as if the particles are "3D"
    double area = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y);
    params.mass_particle = params.rest_density * (area/params.number_fluid_particles_global);

    // Initial spacing between particles
    params.spacing_particle = pow(area/params.number_fluid_particles_global,1.0/2.0);

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
    bytes = params.max_fluid_particles_local * sizeof(fluid_particle);
    total_bytes+=bytes;
    fluid_particle *fluid_particles = malloc(bytes);
    if(fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");
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

    /////////////////////
    // UNIFORM GRID HASH
    /////////////////////
    // +1 added because range begins at 0
    params.grid_size_x = ceil((boundary_global.max_x - boundary_global.min_x) / params.smoothing_radius) + 1;
    params.grid_size_y = ceil((boundary_global.max_y - boundary_global.min_y) / params.smoothing_radius) + 1;
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

    printf("gigabytes allocated: %lld\n", total_bytes/1073741824);

    // Initialize particles
    initParticles(fluid_particle_pointers, fluid_particles, neighbors, hash, &water_volume_global, start_x, number_particles_x, &edges, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, params.number_fluid_particles_local, params.smoothing_radius);

    // Initial configuration
    int fileNum=0;
    writeMPI(fluid_particle_pointers, fileNum++, &params);

    // Main loop
    // In the current form the particles with be re-hashed and halos re-sent for step 0
    int n;
    double start_time, end_time;

    MPI_Barrier(MPI_COMM_WORLD);

    start_time = MPI_Wtime();
    for(n=0; n<params.number_steps; n++) {

//        printf("Rank %d Entering fluid step %d with %d particles\n",rank, n, params.number_fluid_particles_local);

        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

        hash_fluid(fluid_particle_pointers, neighbors, hash, &params);

        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

        hash_halo(fluid_particle_pointers, neighbors, hash, &params);

        // This part is incorrect as halo particles do not have correct pressure/density
        // Need to do a little more mpi work for this to be correct.

        updatePressures(fluid_particle_pointers, neighbors, &params);

        updateAccelerations(fluid_particle_pointers, neighbors, &params);

        updatePositions(fluid_particle_pointers, &out_of_bounds, &edges, &boundary_global, &params);

        if (n % 10 == 0)
            checkPartition(fluid_particle_pointers, &out_of_bounds, &params);

        transferOOBParticles(fluid_particle_pointers, fluid_particles, &out_of_bounds, &params);

        if (n % steps_per_frame == 0)
          writeMPI(fluid_particle_pointers, fileNum++, &params);

    }
    end_time = MPI_Wtime();
    printf("Rank %d Elapsed seconds: %f, num particles: %d\n", rank, end_time-start_time, params.number_fluid_particles_local);

    // Release memory
    free(fluid_particle_pointers);
    free(fluid_particles);
    free(neighbors);
    free(hash);

    // Close MPI
    freeMpiTypes();
    MPI_Finalize();

    return 0;
}

////////////////////////////////////////////////////////////////////////////
// Smoothing Kernels
///////////////////////////////////////////////////////////////////////////

double lap_W_visc(const double r, const double h)
{
    static const double coef = 30.0/M_PI;
    double lap = coef/(h*h*h*h*h) * (h-r);
    return lap;
}

// (h^2 - r^2)^3 normalized in 2D
double W_dens(const double r2, const double h)
{
    const double C = 4.0/(M_PI*h*h*h*h*h*h*h*h);
    double W = C*(h*h-r2)*(h*h-r2)*(h*h-r2);
    return W;

}

// Gradient (h-r)^3 normalized in 2D
double del_W_pressure(const double r, const double h)
{
    static const double coef = -30.0/M_PI;
    double C = coef/(h*h*h*h*h * r);
    double del_W = C*(h-r)*(h-r);
    return  del_W;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////
double computeDensity(const double r2, const param *const params)
{
    double density = params->mass_particle * W_dens(r2,params->smoothing_radius);
    return density;
}

double computePressure(const fluid_particle *const p, const param *const params)
{
    double pressure = 5.0 * (p->density - params->rest_density);
    return pressure;
}

// This should be merged into the neighbor search
void updatePressures(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i, j;
    fluid_particle *p, *q;
    neighbor* n;
    double r2;
    double density;
    double pressure;
    const double h = params->smoothing_radius;

    // This should be moved into the neighbor search
    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        // Self contribution as to not double count below
        p->density = computeDensity(0.0,params);
    }

    int num_neighbors_total = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];
        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
	    r2 = (p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y);
            density = computeDensity(r2,params);
            p->density+=density;
            q->density+=density;
        }
    }

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        p->pressure = computePressure(p,params);
    }

}

// Compute force(accleration) on fluid particle p by fluid particle q
void computeAcceleration(fluid_particle *p, fluid_particle *q, param *params)
{
    const double h = params->smoothing_radius;

    const double x_diff = (p->x - q->x);
    const double y_diff = (p->y - q->y);
    const double mass = params->mass_particle;
    const double q_density = q->density;
    const double p_density = p->density;

    const double r = sqrt(x_diff*x_diff + y_diff*y_diff);

    double accel;
    double a_x, a_y;

    // Pressure force
    if(r>0.0)
        accel = (p->pressure + q->pressure)/(2.0 * q_density)*mass/p_density  * del_W_pressure(r,h);
    // We want some pressure if particles are in same position
    else
        accel = (p->pressure + q->pressure)/(2.0 * q_density)*mass/p_density  * del_W_pressure(h*0.01,h);
    a_x = -accel * x_diff;
    a_y = -accel * y_diff;

    // Viscosity force
    accel = 3.5*mass/q_density/p_density * lap_W_visc(r, h);

    a_x += accel * (q->v_x - p->v_x);
    a_y += accel * (q->v_y - p->v_y);

    // BT 07 http://cg.informatik.uni-freiburg.de/publications/2011_GRAPP_airBubbles.pdf
    //Surface tension
    accel = params->surface_tension * W_dens(r*r,h);
    a_x -= accel * x_diff;
    a_y -= accel * y_diff;

    // Update particle pairs
    p->a_x += a_x;
    p->a_y += a_y;

    q->a_x -= a_x;
    q->a_y -= a_y;

}

// Update particle acclerations
void updateAccelerations(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;

    // This should be reset in neighbor search or somethign
    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];

        p->a_x = 0.0;
        p->a_y = -params->g;
    }

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];

        // Acceleration on p and q due to neighbor fluid particles
        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
            computeAcceleration(p,q,params);
        }
    }
}

// Leap Frog integration with v(t+1) estimated
void updateParticle(fluid_particle *p, int particle_index, param *params)
{
    double dt = params->time_step;

    // Velocity at t + dt/2
    p->v_half_x = p->v_half_x + dt * p->a_x;
    p->v_half_y = p->v_half_y + dt * p->a_y;

    // Velocity at t + dt, must estimate for foce calc
    p->v_x = p->v_half_x + p->a_x * (dt / 2.0);
    p->v_y = p->v_half_y + p->a_y * (dt / 2.0);

    // Position at time t + dt
    p->x = p->x + dt * p->v_half_x;
    p->y = p->y + dt * p->v_half_y;

}

// Update particle position and check boundary
void updatePositions(fluid_particle **fluid_particle_pointers, oob *out_of_bounds, edge *edges, AABB *boundary_global, param *params)
{
    int i;
    fluid_particle *p;
    double h = params->smoothing_radius;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        updateParticle(p, i, params);
        boundaryConditions(p, boundary_global, params);

        // Set OOB particle indicies and update number
        if (p->x < params->node_start_x)
            out_of_bounds->oob_pointer_indicies_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (p->x > params->node_end_x)
            out_of_bounds->oob_pointer_indicies_right[out_of_bounds->number_oob_particles_right++] = i;
    }
}

// Seed simulation with Euler step v(t-dt/2) needed by leap frog integrator
void eulerStart(fluid_particle **fluid_particle_pointers, neighbor* neighbors, param *params)
{
    updatePressures(fluid_particle_pointers, neighbors, params);

    updateAccelerations(fluid_particle_pointers, neighbors, params);

    // Set V (t0 - dt/2)
    int i;
    double dt_half;
    fluid_particle *p;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = fluid_particle_pointers[i];

        dt_half = params->time_step/2.0;
        p->v_half_x = p->v_x - p->a_x * dt_half;
        p->v_half_y = p->v_y - p->a_y * dt_half;
    }
}

// Project particle to AABB surface and reflect velocity
void reflectParticle(fluid_particle *p, param* params, double pen_depth, double *norm)
{
    double dt = params->time_step;
    double dt_half = dt/2.0;
    double restitution = 0.0; // 0 = No slip condition, applicable to viscious fluids

    // project particle back onto surface
    p->x = p->x + (pen_depth * norm[0]);
    p->y = p->y + (pen_depth * norm[1]);

    //////////////////////////////////////
    // Hacky stuff - fix fix fix...please
    if(p->x < 0.0)
	p->x = 0.0;
    else if(p->x > 1.0)
	p->x = 1.0;
    if(p->y < 0.0)
	p->y = 0.0;
    else if(p->y > 1.0)
	p->y = 1.0;

     norm[0] = (double)sgn(norm[0]);
     norm[1] = (double)sgn(norm[1]);
     double mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
    norm[0]/=mag;
    norm[1]/=mag;
    ///////////////////////////////////////


    double v_mag = sqrt(p->v_x*p->v_x + p->v_y*p->v_y);
    double vDotn = p->v_x*norm[0] + p->v_y*norm[1];

    double v_half_mag = sqrt(p->v_half_x*p->v_half_x + p->v_half_y*p->v_half_y);
    double vHalfDotn = p->v_half_x*norm[0] + p->v_half_y*norm[1];

    // Calculate new reflected velocity
    p->v_x = p->v_x - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[0]);
    p->v_y = p->v_y - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[1]);

    p-> v_half_x = p-> v_half_x - ((1.0 + (restitution*pen_depth)/(dt/2.0*v_half_mag)) * vHalfDotn * norm[0]);
    p-> v_half_y = p-> v_half_y - ((1.0 + (restitution*pen_depth)/(dt/2.0*v_half_mag)) * vHalfDotn * norm[1]);
}

// Assume AABB with min point being axis origin
void boundaryConditions(fluid_particle *p, AABB *boundary, param *params)
{
    // AABB center
    double c_x = boundary->min_x + (boundary->max_x - boundary->min_x)/2.0;
    double c_y = boundary->min_y + (boundary->max_y - boundary->min_y)/2.0;

    // AABB extent from center of frame, aka half size
    double e_x = boundary->max_x - c_x;
    double e_y = boundary->max_y - c_y;

    // local position of particle from center of AABB
    double local_x = p->x - c_x;
    double local_y = p->y - c_y;

    // Relative distance between particle and boundary
    double f_x = fabs(local_x) - e_x;
    double f_y = fabs(local_y) - e_y;

    // F = max(f_x,f_y,f_z): F > 0 outside, F < 0 inside, F=0 on a face
    double f = max(f_x,f_y);

    // Particle outside of bounding volume and should be inside
    if(f > 0.0)
    {
        // Calculate local contact point
        double local_cp_x = min(e_x, max(-e_x, local_x));
        double local_cp_y = min(e_y, max(-e_y, local_y));

        // world position of contact point // Must modify for OBB
        double cp_x = c_x + local_cp_x;
        double cp_y = c_y + local_cp_y;

        // Penetration depth
        double depth = sqrt((cp_x-p->x)*(cp_x-p->x) + (cp_y-p->y)*(cp_y-p->y));

        // Surface normal scaled
        // Use local coordinates otherwise small round off errors from world->local frame throw off sgn calculation
//        double n_x = (double)sgn(local_cp_x - local_x);
//        double n_y = (double)sgn(local_cp_y - local_y);
//        double n_z = (double)sgn(local_cp_z - local_z);
        double n_x = (local_cp_x - local_x);
        double n_y = (local_cp_y - local_y);


        double norm = sqrt(n_x*n_x + n_y*n_y);
        n_x = n_x/norm;
        n_y = n_y/norm;

        double normal[2] = {n_x,n_y};
        reflectParticle(p, params, depth, normal);
    }
    else {

    //////////////////////////////////////
    // Hacky stuff - fix fix fix...please
    if(p->x < 0.0)
        p->x = 0.0;
    else if(p->x > 1.0)
        p->x = 1.0;
    if(p->y < 0.0)
        p->y = 0.0;
    else if(p->y > 1.0)
        p->y = 1.0;


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
        fluid_particle_pointers[i]->density = params->rest_density;
    }

    // Send halo particles
    startHaloExchange(fluid_particle_pointers,fluid_particles, edges, params);

    //Generate neighbor hash
    hash_fluid(fluid_particle_pointers, neighbors, hash, params);
    finishHaloExchange(fluid_particle_pointers,fluid_particles, edges, params);
    hash_halo(fluid_particle_pointers, neighbors, hash, params);

    // Set initial velocity half step
    eulerStart(fluid_particle_pointers, neighbors, params);
}
