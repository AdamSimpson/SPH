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
    params.number_steps = 1000;
    params.time_step = 1/60.0;
    params.number_fluid_particles_global = 1800;
    params.rest_density = 1000.0; // water: kg/m^3
    
    // Boundary box
    boundary_global.min_x = 0.0;
    boundary_global.max_x = 0.5;
    boundary_global.min_y = 0.0;
    boundary_global.max_y = 0.5;
    boundary_global.min_z = 0.0;
    boundary_global.max_z = 0.5;
    
    // water volume
    water_volume_global.min_x = 0.005;
    water_volume_global.max_x = 0.495;
    water_volume_global.min_y = 0.005;
    water_volume_global.max_y = 0.495;
    water_volume_global.min_z = 0.005;
    water_volume_global.max_z = 0.4;
    
    // Mass of each particle
    double volume = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y) * (water_volume_global.max_z - water_volume_global.min_z);
    params.mass_particle = params.rest_density * (volume/params.number_fluid_particles_global);
    
    // Initial spacing between particles
    params.spacing_particle = pow(volume/params.number_fluid_particles_global,1.0/3.0);

    // Smoothing radius, h
    params.smoothing_radius = 2.0*params.spacing_particle;

    printf("smoothing radius: %f\n", params.smoothing_radius);    

    // Number of steps before frame needs to be written for 30 fps
    int steps_per_frame = (int)(1.0/(params.time_step*30.0));
    
    // Calculate speed of sound for simulation
    double max_height = water_volume_global.max_y;
    double max_velocity = sqrt(2.0*params.g*max_height);
    params.speed_sound = max_velocity/sqrt(0.01);

    double recomend_step = 0.4 * params.smoothing_radius / (params.speed_sound * (1+ 0.6*params.alpha));
    printf("Using time step: %f, Minimum recomended %f\n",params.time_step, recomend_step);    

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

    // Allocate new hash with all values zeroed
    params.length_hash = 2 * (params.max_fluid_particles_local); // Hash works best with nearest prime > 2*nt
    n_bucket* hash = calloc(params.length_hash, sizeof(n_bucket));
    total_bytes+=params.length_hash*sizeof(n_bucket);
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
    free(fluid_particles);
    free(fluid_particle_pointers);
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
double lap_W_visc(fluid_particle *p, fluid_particle *q, double r, double h)
{
    static const coef = 45.0/M_PI;
    double lap = coef/(h*h*h*h*h*h) * (h-r);
    return lap;
}

double W_dens(fluid_particle *p, fluid_particle *q, double r, double h)
{
    const static coef = 315.0/(64.0 * 3.14);
    double C = coef/(h*h*h*h*h*h*h*h*h);
    double W = C*(h*h-r*r)*(h*h-r*r)*(h*h-r*r);
    return W;

}

double del_W_pressure(fluid_particle *p, fluid_particle *q, double r, double h)
{
    const static coef = -45.0/M_PI;
    double C = coef/(h*h*h*h*h*h * r);
    double del_W = C*(h-r)*(h-r);
    return  del_W;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////
double computeDensity(fluid_particle *p, fluid_particle *q, double r, param *params)
{
    double density = params->mass_particle * W_dens(p,q,r,params->smoothing_radius);
    return density;
}

double computePressure(fluid_particle *p, param *params)
{
    double pressure = 3.0 * (p->density - params->rest_density);
    return pressure;
}

void updatePressures(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i, j;
    fluid_particle *p, *q;
    neighbor* n;
    double r;
    const double h = params->smoothing_radius;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];
 	p->density = 0.0;       
        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
            r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y) + (p->z-q->z)*(p->z-q->z));
	    if(r <= h){
                q = n->fluid_neighbors[j];
                p->density += computeDensity(p,q,r,params);
	    }
        }
        p->pressure = computePressure(p,params);
    }
}

// Compute force(accleration) on fluid particle p by fluid particle q
void computeAcceleration(fluid_particle *p, fluid_particle *q, param *params)
{
    double accel;
    const double h = params->smoothing_radius;
    const double x_diff = (p->x - q->x);
    const double y_diff = (p->y - q->y);
    const double z_diff = (p->z - q->z);
    const double r = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
    
    if(p == q || r > h)
	return;

    // Pressure force
    accel = (p->pressure + q->pressure)/(2.0 * q->density)*params->mass_particle/p->density  * del_W_pressure(p,q,r,h);
    p->a_x -= accel * x_diff;
    p->a_y -= accel * y_diff;
    p->a_z -= accel * z_diff;
        
    // Viscosity force
    double mu = 3.5;
    accel = mu*params->mass_particle/q->density/p->density * lap_W_visc(p, q, r, h);
    p->a_x += accel  * (q->v_x - p->v_x);
    p->a_y += accel  * (q->v_y - p->v_y);
    p->a_z += accel  * (q->v_z - p->v_z);

    // BT 07 http://cg.informatik.uni-freiburg.de/publications/2011_GRAPP_airBubbles.pdf
    //Surface tension
    accel = params->surface_tension * W_dens(p,q,r,h);
    p->a_x -= accel * x_diff;
    p->a_y -= accel * y_diff;
    p->a_z -= accel * z_diff;

}

// Update particle acclerations
void updateAccelerations(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;
    
    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        
        n = &neighbors[i];
        
        p->a_x = 0.0;
        p->a_y = 0.0;
        p->a_z = -params->g;
        
        // Acceleration on p due to neighbor fluid particles
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
    p->v_half_z = p->v_half_z + dt * p->a_z;
    
    // Velocity at t + dt, must estimate for foce calc
    p->v_x = p->v_half_x + p->a_x * (dt / 2.0);
    p->v_y = p->v_half_y + p->a_y * (dt / 2.0);
    p->v_z = p->v_half_z + p->a_z * (dt / 2.0);
    
    // Position at time t + dt
    p->x = p->x + dt * p->v_half_x;
    p->y = p->y + dt * p->v_half_y;
    p->z = p->z + dt * p->v_half_z;

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
        p->v_half_z = p->v_z - p->a_z * dt_half;
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
    p->z = p->z + (pen_depth * norm[2]);

    //////////////////////////////////////
    // Hacky stuff - fix fix fix...please
    if(p->x < 0.0)
	p->x = 0.0;

     norm[0] = (double)sgn(norm[0]);
     norm[1] = (double)sgn(norm[1]);
     norm[2] = (double)sgn(norm[2]);
     double mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1] + norm[2]*norm[2]);
    norm[0]/=mag;
    norm[1]/=mag;
    norm[2]/=mag;
    ///////////////////////////////////////


    double v_mag = sqrt(p->v_x*p->v_x + p->v_y*p->v_y + p->v_z*p->v_z);
    double vDotn = p->v_x*norm[0] + p->v_y*norm[1] + p->v_z*norm[2];

    double v_half_mag = sqrt(p->v_half_x*p->v_half_x + p->v_half_y*p->v_half_y + p->v_half_z*p->v_half_z);
    double vHalfDotn = p->v_half_x*norm[0] + p->v_half_y*norm[1] + p->v_half_z*norm[2];

    // Calculate new reflected velocity
    p->v_x = p->v_x - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[0]);
    p->v_y = p->v_y - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[1]);
    p->v_z = p->v_z - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[2]);

    p-> v_half_x = p-> v_half_x - ((1.0 + (restitution*pen_depth)/(dt/2.0*v_half_mag)) * vHalfDotn * norm[0]);
    p-> v_half_y = p-> v_half_y - ((1.0 + (restitution*pen_depth)/(dt/2.0*v_half_mag)) * vHalfDotn * norm[1]);
    p-> v_half_z = p-> v_half_z - ((1.0 + (restitution*pen_depth)/(dt/2.0*v_half_mag)) * vHalfDotn * norm[2]);
}

// Assume AABB with min point being axis origin 
void boundaryConditions(fluid_particle *p, AABB *boundary, param *params)
{
    // AABB center
    double c_x = boundary->min_x + (boundary->max_x - boundary->min_x)/2.0;
    double c_y = boundary->min_y + (boundary->max_y - boundary->min_y)/2.0;
    double c_z = boundary->min_z + (boundary->max_z - boundary->min_z)/2.0;

    // AABB extent from center of frame, aka half size
    double e_x = boundary->max_x - c_x;
    double e_y = boundary->max_y - c_y;
    double e_z = boundary->max_z - c_z;

    // local position of particle from center of AABB
    double local_x = p->x - c_x;
    double local_y = p->y - c_y;
    double local_z = p->z - c_z;

    // Relative distance between particle and boundary
    double f_x = fabs(local_x) - e_x;
    double f_y = fabs(local_y) - e_y;
    double f_z = fabs(local_z) - e_z;

    // F = max(f_x,f_y,f_z): F > 0 outside, F < 0 inside, F=0 on a face
    double f = max(max(f_x,f_y),f_z);

    // Particle outside of bounding volume and should be inside
    if(f > 0.0)
    {
        // Calculate local contact point
        double local_cp_x = min(e_x, max(-e_x, local_x));
        double local_cp_y = min(e_y, max(-e_y, local_y));
        double local_cp_z = min(e_z, max(-e_z, local_z));

        // world position of contact point // Must modify for OBB
        double cp_x = c_x + local_cp_x;
        double cp_y = c_y + local_cp_y;
        double cp_z = c_z + local_cp_z;

        // Penetration depth
        double depth = sqrt((cp_x-p->x)*(cp_x-p->x) + (cp_y-p->y)*(cp_y-p->y) + (cp_z-p->z)*(cp_z-p->z));

        // Surface normal scaled
        // Use local coordinates otherwise small round off errors from world->local frame throw off sgn calculation
//        double n_x = (double)sgn(local_cp_x - local_x);
//        double n_y = (double)sgn(local_cp_y - local_y);
//        double n_z = (double)sgn(local_cp_z - local_z);
        double n_x = (local_cp_x - local_x);
        double n_y = (local_cp_y - local_y);
        double n_z = (local_cp_z - local_z);


        double norm = sqrt(n_x*n_x + n_y*n_y + n_z*n_z);
        n_x = n_x/norm;
        n_y = n_y/norm;
        n_z = n_z/norm;

        double normal[3] = {n_x,n_y,n_z};
        reflectParticle(p, params, depth, normal);
    }
}

// Initialize particles
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                neighbor *neighbors, n_bucket *hash, AABB* water, int start_x, int number_particles_x, edge *edges, param* params)
{
    int i;
    fluid_particle *p;
    double h = params->spacing_particle;

    // Create fluid volume
    constructFluidVolume(fluid_particle_pointers, fluid_particles, water, start_x, number_particles_x, edges, params);

    // NULL out unused fluid pointers
    for(i=params->number_fluid_particles_local; i<params->max_fluid_particles_local; i++)
        fluid_particle_pointers[i] = NULL;

    // Initialize particle values
    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particle_pointers[i]->a_x = 0.0;
        fluid_particle_pointers[i]->a_y = 0.0;
        fluid_particle_pointers[i]->a_z = 0.0;
        fluid_particle_pointers[i]->v_x = 0.0;
        fluid_particle_pointers[i]->v_y = 0.0;
        fluid_particle_pointers[i]->v_z = 0.0;
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
