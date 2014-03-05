#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

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

    printf("Global rank %d\n", rank);
    MPI_Finalize();
    return 0;
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

    unsigned int i;

    params.tunable_params.kill_sim = false;
    params.tunable_params.g = 6.0f;
    params.tunable_params.time_step = 1.0f/30.0f;
    params.tunable_params.k = 0.2f;
    params.tunable_params.k_near = 6.0f;
    params.tunable_params.k_spring = 10.0f;
    params.tunable_params.sigma = 20.0f;
    params.tunable_params.beta = 2.0f;
    params.tunable_params.rest_density = 30.0f;
    params.tunable_params.mover_width = 1.0f;
    params.tunable_params.mover_height = 1.0f;
    params.tunable_params.mover_type = SPHERE_MOVER;

    // The number of particles used may differ slightly
    params.number_fluid_particles_global = 1000;

    // Boundary box
    // This simulation assumes in various spots min is 0.0
    boundary_global.min_x = 0.0f;
    boundary_global.max_x = 15.0f;
    boundary_global.min_y = 0.0f;

    // Receive aspect ratio to scale world y max
    short pixel_dims[2];
    float aspect_ratio;
    MPI_Bcast(pixel_dims, 2, MPI_SHORT, 0, MPI_COMM_WORLD);
    aspect_ratio = (float)pixel_dims[0]/(float)pixel_dims[1];
    boundary_global.max_y = boundary_global.max_x / aspect_ratio;
    
    // water volume
    water_volume_global.min_x = 0.0f;
    water_volume_global.max_x = boundary_global.max_x;
    water_volume_global.min_y = 0.0f;
    water_volume_global.max_y = boundary_global.max_y;

    params.number_halo_particles = 0;

    int start_x;  // where in x direction this nodes particles start
    int number_particles_x; // number of particles in x direction for this node

    // Fluid area in initial configuration
    float area = (water_volume_global.max_x - water_volume_global.min_x) * (water_volume_global.max_y - water_volume_global.min_y);

    // Initial spacing between particles
    float spacing_particle = pow(area/params.number_fluid_particles_global,1.0/2.0);

    // Divide problem set amongst nodes
    partitionProblem(&boundary_global, &water_volume_global, &start_x, &number_particles_x, spacing_particle, &params);

    // Set local/global number of particles to allocate
    setParticleNumbers(&boundary_global, &water_volume_global, &edges, &out_of_bounds, number_particles_x, spacing_particle, &params);

    // We will allocate enough room for all particles on single node
    // We also must take into account halo particles are placed onto the end of the max particle index
    // So this value can be even greater than the number of global
    // Before reaching this point the program should, but doesn't, intelligenly clean up fluid_particles
    int max_fluid_particles_local = 2*params.number_fluid_particles_global;

    // Smoothing radius, h
    params.tunable_params.smoothing_radius = 2.0f*spacing_particle;

    printf("smoothing radius: %f\n", params.tunable_params.smoothing_radius);

    // Send initial world dimensions and max particle count to render node
    if(rank == 0) {
        float world_dims[2];
        world_dims[0] = boundary_global.max_x;
        world_dims[1] = boundary_global.max_y;
        MPI_Send(world_dims, 2, MPI_FLOAT, 0, 8, MPI_COMM_WORLD);
	MPI_Send(&params.number_fluid_particles_global, 1, MPI_INT, 0, 9, MPI_COMM_WORLD);
    }

    // Neighbor grid setup
    neighbor_grid_t neighbor_grid;
    neighbor_grid.max_bucket_size = 100;
    neighbor_grid.max_neighbors = neighbor_grid.max_bucket_size*4;
    neighbor_grid.spacing = params.tunable_params.smoothing_radius;

    size_t total_bytes = 0;
    size_t bytes;
    // Allocate fluid particles array
    bytes = max_fluid_particles_local * sizeof(fluid_particle);
    total_bytes+=bytes;
    fluid_particle *fluid_particles = malloc(bytes);
    if(fluid_particles == NULL)
        printf("Could not allocate fluid_particles\n");

    // Allocate (x,y) coordinate array, transfer pixel coords
    bytes = 2 * max_fluid_particles_local * sizeof(short);
    total_bytes+=bytes;
    short *fluid_particle_coords = malloc(bytes);
    if(fluid_particle_coords == NULL)
        printf("Could not allocate fluid_particle coords\n");

    // Allocate pointer array used to traverse non vacant particles
    bytes = max_fluid_particles_local * sizeof(fluid_particle*);
    total_bytes+=bytes;
    fluid_particle **fluid_particle_pointers = malloc(bytes);
    if(fluid_particle_pointers == NULL)
        printf("Could not allocate fluid_particle_pointers\n");

    // Allocate neighbor array
    neighbor *neighbors = calloc(max_fluid_particles_local, sizeof(neighbor));
    fluid_particle **fluid_neighbors = calloc(max_fluid_particles_local * neighbor_grid.max_neighbors, sizeof(fluid_particle *));
    // Set pointer in each bucket
    for(i=0; i< max_fluid_particles_local; i++ )
        neighbors[i].fluid_neighbors = &(fluid_neighbors[i*neighbor_grid.max_neighbors]);

    neighbor_grid.neighbors = neighbors;
    total_bytes+= (max_fluid_particles_local*sizeof(neighbor) + neighbor_grid.max_neighbors*sizeof(fluid_particle *));
    if(neighbors == NULL || fluid_neighbors == NULL)
        printf("Could not allocate neighbors\n");

    // UNIFORM GRID HASH
    neighbor_grid.size_x = ceil((boundary_global.max_x - boundary_global.min_x) / neighbor_grid.spacing);
    neighbor_grid.size_y = ceil((boundary_global.max_y - boundary_global.min_y) / neighbor_grid.spacing);
    unsigned int length_hash = neighbor_grid.size_x * neighbor_grid.size_y;
    printf("grid x: %d grid y %d\n", neighbor_grid.size_x, neighbor_grid.size_y);
    bucket_t* grid_buckets = calloc(length_hash, sizeof(bucket_t));
    fluid_particle **bucket_particles = calloc(length_hash * neighbor_grid.max_bucket_size, sizeof(fluid_particle *));
    neighbor_grid.grid_buckets = grid_buckets;
    for(i=0; i < length_hash; i++)
	grid_buckets[i].fluid_particles = &(bucket_particles[i*neighbor_grid.max_bucket_size]);
    total_bytes+= (length_hash * sizeof(bucket_t) + neighbor_grid.max_bucket_size * sizeof(fluid_particle *));
    if(grid_buckets == NULL || bucket_particles == NULL)
        printf("Could not allocate hash\n");

    // Allocate edge index arrays
    edges.edge_pointers_left = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    edges.edge_pointers_right = malloc(edges.max_edge_particles * sizeof(fluid_particle*));
    // Allocate out of bound index arrays
    out_of_bounds.oob_pointer_indicies_left = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.oob_pointer_indicies_right = malloc(out_of_bounds.max_oob_particles * sizeof(int));
    out_of_bounds.vacant_indicies = malloc(2*out_of_bounds.max_oob_particles * sizeof(int));

    printf("bytes allocated: %lu\n", total_bytes);

    // Initialize particles
    initParticles(fluid_particle_pointers, fluid_particles, &water_volume_global, start_x,
		  number_particles_x, &edges, max_fluid_particles_local, spacing_particle, &params);

    // Print some parameters
    printf("Rank: %d, fluid_particles: %d, smoothing radius: %f \n", rank, params.number_fluid_particles_local, params.tunable_params.smoothing_radius);

    // Send intiial paramaters to render node
    tunable_parameters *null_tunable_param = NULL;
    int *null_recvcnts = NULL;
    int *null_displs = NULL;
    MPI_Gatherv(&params.tunable_params, 1, TunableParamtype, null_tunable_param, null_recvcnts, null_displs, TunableParamtype, 0, MPI_COMM_WORLD);

    fluid_particle *p;
    unsigned int n = 0;
    fluid_particle *null_particle = NULL;
    float *null_float = NULL;

    MPI_Request coords_req = MPI_REQUEST_NULL;

    // Main simulation loop
    while(1) {

        // Initialize velocities
        apply_gravity(fluid_particle_pointers, &params);

        // Viscosity impluse
        viscosity_impluses(fluid_particle_pointers, neighbors, &params);

        // Advance to predicted position and set OOB particles
        predict_positions(fluid_particle_pointers, &boundary_global, &params);

	    // Make sure that async send to render node is complete
        if(coords_req != MPI_REQUEST_NULL)
	    MPI_Wait(&coords_req, MPI_STATUS_IGNORE);

        // Receive updated paramaters from render nodes
        MPI_Scatterv(null_tunable_param, 0, null_displs, TunableParamtype, &params.tunable_params, 1, TunableParamtype, 0,  MPI_COMM_WORLD);

	if(params.tunable_params.kill_sim)
	    break;

        // Identify out of bounds particles and send them to appropriate rank
        identify_oob_particles(fluid_particle_pointers, fluid_particles, &out_of_bounds, &boundary_global, &params);

        // Hash the non halo regions
        // This will update the densities so when the halo is exchanged the halo particles are up to date
        // This works well on the raspi's but destroys communication/computation overlap
        hash_fluid(fluid_particle_pointers, &neighbor_grid, &params, true);

         // Exchange halo particles
        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);
        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

	    // Add the halo particles to neighbor buckets
	    // Also update density
        hash_halo(fluid_particle_pointers, &neighbor_grid, &params, true);

        // double density relaxation
	    // halo particles will be missing origin contributions to density/pressure
        double_density_relaxation(fluid_particle_pointers, neighbors, &params);

        // update velocity
        updateVelocities(fluid_particle_pointers, &edges, &boundary_global, &params);

        // Not updating halo particles and hash after relax can be used to speed things up
        // Not updating these can cause unstable behavior

        // Exchange halo particles from relaxed positions
//        startHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

        // We can hash during exchange as the density is not needed
        hash_fluid(fluid_particle_pointers, &neighbor_grid, &params, false);

        // Finish asynch halo exchange
//        finishHaloExchange(fluid_particle_pointers,fluid_particles, &edges, &params);

        // Update hash with relaxed positions
//        hash_halo(fluid_particle_pointers, &neighbor_grid, &params, false);

        // We do not transfer particles that have gone OOB since relaxation
        // to reduce communication cost

        // Pack fluid particle coordinates
	// This should be sent as short in pixel coordinates
        for(i=0; i<params.number_fluid_particles_local; i++) {
            p = fluid_particle_pointers[i];
            fluid_particle_coords[i*2] = (2.0f*p->x/boundary_global.max_x - 1.0f) * SHRT_MAX; // convert to short using full range
            fluid_particle_coords[(i*2)+1] = (2.0f*p->y/boundary_global.max_y - 1.0f) * SHRT_MAX; // convert to short using full range
        }
	// Async send fluid particle coordinates to render node
	MPI_Isend(fluid_particle_coords, 2*params.number_fluid_particles_local, MPI_SHORT, 0, 17, MPI_COMM_WORLD, &coords_req);

        // iterate sim loop counter
        n++;
    }

    // Release memory
/*
    free(fluid_particle_pointers);
    free(fluid_particles);
    free(neighbors);
    free(hash);

    // Close MPI
    freeMpiTypes();
*/
}

// This should go into the hash, perhaps with the viscocity?
void apply_gravity(fluid_particle **fluid_particle_pointers, param *params)
{
    int i;
    fluid_particle *p;
    float dt = params->tunable_params.time_step;
    float g = -params->tunable_params.g;

    for(i=0; i<(params->number_fluid_particles_local + params->number_halo_particles); i++) {
        p = fluid_particle_pointers[i];
        p->v_y += g*dt;

        // Zero out density as well
        p->density = 0.0f;
        p->density_near = 0.0f;
     }
}

// Add viscosity impluses
void viscosity_impluses(fluid_particle **fluid_particle_pointers, neighbor* neighbors, param *params)
{
    int i, j, num_fluid;
    fluid_particle *p, *q;
    neighbor* n;
    float r, r_recip, ratio, u, imp, imp_x, imp_y;
    float p_x, p_y;
    float QmP_x, QmP_y;
    float h_recip, sigma, beta, dt;

    num_fluid = params->number_fluid_particles_local;
    h_recip = 1.0f/params->tunable_params.smoothing_radius;
    sigma = params->tunable_params.sigma;
    beta = params->tunable_params.beta;
    dt = params->tunable_params.time_step;


    for(i=params->number_fluid_particles_local; i-- > 0; ) {
//    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];
 	    p_x = p->x;
	    p_y = p->y;

        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
	
            QmP_x = (q->x-p_x);
            QmP_y = (q->y-p_y);
            r = sqrt(QmP_x*QmP_x + QmP_y*QmP_y);

            r_recip = 1.0f/r;
            ratio = r*h_recip;

            //Inward radial velocity
            u = ((p->v_x-q->v_x)*QmP_x + (p->v_y-q->v_y)*QmP_y)*r_recip;
            if(u>0.0f)
            {
                imp = dt * (1-ratio)*(sigma * u + beta * u*u);
                imp_x = imp*QmP_x*r_recip;
                imp_y = imp*QmP_y*r_recip;

		// Not correct to use velocity check but will stop velocity from
		// blowing up
		checkVelocity(&imp_x, &imp_y);

                p->v_x -= imp_x*0.5f;
                p->v_y -= imp_y*0.5f;

                if(q->id < num_fluid) {
                    q->v_x += imp_x*0.5f;
                    q->v_y += imp_y*0.5f;

                }
                else { // Only apply half of the impulse to halo particles as they are missing "home" contribution
                    q->v_x += imp_x*0.125f;
                    q->v_y += imp_y*0.125f;
                }
                
            }

        }
    }
}

// Identify out of bounds particles and send them to appropriate rank
void identify_oob_particles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles, oob *out_of_bounds, AABB *boundary_global, param *params)
{
    int i;
    fluid_particle *p;

    // Reset OOB numbers
    out_of_bounds->number_oob_particles_left = 0;
    out_of_bounds->number_oob_particles_right = 0;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];

        // Set OOB particle indicies and update number
        if (p->x < params->tunable_params.node_start_x)
            out_of_bounds->oob_pointer_indicies_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (p->x >= params->tunable_params.node_end_x)
            out_of_bounds->oob_pointer_indicies_right[out_of_bounds->number_oob_particles_right++] = i;
    }
 
   // Transfer particles that have left the processor bounds
   transferOOBParticles(fluid_particle_pointers, fluid_particles, out_of_bounds, params);
}



// Predict position
void predict_positions(fluid_particle **fluid_particle_pointers, AABB *boundary_global, param *params)
{
    int i;
    fluid_particle *p;
    float dt = params->tunable_params.time_step;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        p = fluid_particle_pointers[i];
	p->x_prev = p->x;
        p->y_prev = p->y;
	p->x += (p->v_x * dt);
        p->y += (p->v_y * dt);

	// Enforce boundary conditions
        boundaryConditions(p, boundary_global, params);
    }
}

// Calculate the density contribution of p on q and q on p
// r is passed in as this function is called in the hash which must also calculate r
void calculate_density(fluid_particle *p, fluid_particle *q, float ratio)
{

    float OmR2 = (1.0f-ratio)*(1.0f-ratio); // (one - r)^2
    if(ratio < 1.0f) {
	p->density += OmR2;
	p->density_near += OmR2*(1.0f-ratio);

	q->density += OmR2;
	q->density_near += OmR2*(1.0f-ratio);
    }

}

void double_density_relaxation(fluid_particle **fluid_particle_pointers, neighbor *neighbors, param *params)
{
    int i, j, num_fluid;
    fluid_particle *p, *q;
    neighbor* n;
    float r,ratio,dt,h,h_recip,r_recip,D,D_x,D_y;
    float k, k_near, k_spring, p_pressure, p_pressure_near, rest_density;
    float OmR;

    num_fluid = params->number_fluid_particles_local;
    k = params->tunable_params.k;
    k_near = params->tunable_params.k_near;
    k_spring = params->tunable_params.k_spring;
    h = params->tunable_params.smoothing_radius;
    h_recip = 1.0f/h;
    dt = params->tunable_params.time_step;
    rest_density = params->tunable_params.rest_density;

    // Calculate the pressure of all particles, including halo
    for(i=0; i<params->number_fluid_particles_local + params->number_halo_particles; i++) {
        p = fluid_particle_pointers[i];
        // Compute pressure and near pressure
        p->pressure = k * (p->density - rest_density);
        p->pressure_near = k_near * p->density_near;
    }

    // Iterating through the array in reverse reduces biased particle movement
    for(i=params->number_fluid_particles_local; i-- > 0; ) {
        p = fluid_particle_pointers[i];
        n = &neighbors[i];
        p_pressure = p->pressure;
        p_pressure_near = p->pressure_near;

        for(j=0; j<n->number_fluid_neighbors; j++) {

            q = n->fluid_neighbors[j];
            r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y));
	        r_recip = 1.0f/r;
	        ratio = r*h_recip;
	        OmR = 1.0f - ratio;

            // Attempt to move clustered particles apart
            if(r <= 0.000001f) {
                p->x += 0.000001f;
                p->y += 0.000001f;
            }

	    if(ratio < 1.0f && r > 0.0f) {
                // Updating both neighbor pairs at the same time, slightly different than the paper but quicker
                // Also the running sum of D for particle p seems to produce more bias/instability so is removed
                D = dt*dt*((p_pressure+q->pressure)*OmR + (p_pressure_near+q->pressure_near)*OmR*OmR + k_spring*(h-r)*0.5);
                D_x = D*(q->x-p->x)*r_recip;
                D_y = D*(q->y-p->y)*r_recip;

                // Do not move the halo particles full D
                // Halo particles are missing D from their origin so I believe this is appropriate
                if(q->id < num_fluid) {
                    q->x += D_x;
                    q->y += D_y;
                }	
                else { // Move the halo particles only half way to account for other sides missing contribution
                    q->x += D_x*0.125f;
                    q->y += D_y*0.125f;
                }
 
                p->x -= D_x;
                p->y -= D_y;
           }
       }
    }
}

void checkVelocity(float *v_x, float *v_y)
{
    const float v_max = 5.0f;

    if(*v_x > v_max)
        *v_x = v_max;
    else if(*v_x < -v_max)
        *v_x = -v_max;
    if(*v_y > v_max)
        *v_y = v_max;
    else if(*v_y < -v_max)
        *v_y = -v_max;
}

void updateVelocity(fluid_particle *p, param *params)
{
    float dt = params->tunable_params.time_step;
    float v_x, v_y;

    v_x = (p->x-p->x_prev)/dt;
    v_y = (p->y-p->y_prev)/dt;

    checkVelocity(&v_x, &v_y);

    p->v_x = v_x;
    p->v_y = v_y;
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

void collisionImpulse(fluid_particle *p, float norm_x, float norm_y, param *params)
{
    // Boundary friction
    const static float mu = 1.0f;
    float dt = params->tunable_params.time_step;

    float v_dot_n, vx_norm, vy_norm, vx_tan, vy_tan, I_x, I_y;

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

//    printf("pos(%f,%f) norm(%f,%f) V(%f,%f) I(%f,%f)\n", p->x, p->y, norm_x, norm_y, p->v_x, p->v_y,I_x, I_y);

    // Modify particle position

    if(signbit(p->v_x) ==  signbit(norm_x))
        I_x = -I_x;
    if(signbit(p->v_y) ==  signbit(norm_y))
        I_y = -I_y;

    p->x_prev = p->x;
    p->y_prev = p->y;

/*
    p->v_x -= I_x;
    p->v_y -= I_y;

    p->x += p->v_x * dt;
    p->y += p->v_y * dt;
*/

    p->x -= I_x * dt;
    p->y -= I_y * dt;

//   printf("after pos(%f,%f)\n", p->x, p->y);
}

// Assume AABB with min point being axis origin
void boundaryConditions(fluid_particle *p, AABB *boundary, param *params)
{

    float center_x = params->tunable_params.mover_center_x;
    float center_y = params->tunable_params.mover_center_y;

    // Boundary condition for sphere mover
    if(params->tunable_params.mover_type == SPHERE_MOVER)
    {
        // Sphere width == height
        float radius = params->tunable_params.mover_width*0.5f;
        float norm_x;
        float norm_y;

        // Both circle tests can be combined if no impulse is used
        // Test if inside of circle
        float d;
        float d2 = (p->x - center_x)*(p->x - center_x) + (p->y - center_y)*(p->y - center_y);
        if(d2 <= radius*radius && d2 > 0.0f) {
            d = sqrt(d2);
            norm_x = (center_x-p->x)/d;
            norm_y = (center_y-p->y)/d;
	    
	    // With no collision impulse we can handle penetration here
            float pen_dist = radius - d;
            p->x -= pen_dist * norm_x;
            p->y -= pen_dist * norm_y;

           // collisionImpulse(p, norm_x, norm_y, params);
        }

        // Make sure particle is outside of circle
        /*
        d2 = (p->x - center_x)*(p->x - center_x) + (p->y - center_y)*(p->y - center_y);
        if(d2 <= radius*radius) {
            float pen_dist = radius - sqrt(d2);
            p->x -= pen_dist * norm_x;
            p->y -= pen_dist * norm_y;
        }
        */
    }

    // Boundary condition for rectangle mover
    if(params->tunable_params.mover_type == RECTANGLE_MOVER)
    {
        float half_width = params->tunable_params.mover_width*0.5;
        float half_height = params->tunable_params.mover_height*0.5;

        // Particle possition relative to mover center
        float pos_center_x = p->x - center_x;
        float pos_center_y = p->y - center_y;

        // Distance from particle to mover center
	float dist_center_x = fabs(pos_center_x);
	float dist_center_y = fabs(pos_center_y);  

	// Test if inside rectangle
        if( dist_center_x < half_width && dist_center_y < half_height)
        {
            // To find where penetrated from we assume
            // particle is closest to penetrated side

            // Particle penetration depth into rectangle
            float pen_depth_x = half_width - dist_center_x;
            float pen_depth_y = half_height - dist_center_y;

            // Particle closer to left/right sides
            if(pen_depth_x < pen_depth_y){
                // Entered left side
                if(pos_center_x < 0.0f)
                    p->x -= pen_depth_x;
                else // Entered right side
                    p->x += pen_depth_x;
            }
            else { // Particle closer to top/bottom
                // Entered bottom
                if(pos_center_y < 0.0f)
                    p->y -= pen_depth_y;
                else // Entered top
                    p->y += pen_depth_y;
            }
        }
    }
 
    // Boundary seems more stable without collision impulse
    /*
    // Update velocity
    if(p->x <= boundary->min_x) {
        collisionImpulse(p,1.0,0.0,params);
    }
    else if(p->x >= boundary->max_x){
        collisionImpulse(p, -1.0, 0.0, params);
    }
    if(p->y <= boundary->min_y) {
        collisionImpulse(p,0.0,1.0,params);
    }
    else if(p->y >= boundary->max_y){
        collisionImpulse(p,0.0,-1.0,params);
    }
    */

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(p->x < boundary->min_x) {
        p->x = boundary->min_x;
    }
    else if(p->x > boundary->max_x){
        p->x = boundary->max_x-0.001f;
    }
    if(p->y < boundary->min_y) {
        p->y = boundary->min_y;
    }
    else if(p->y > boundary->max_y){
        p->y = boundary->max_y-0.001f;
    }

}

// Initialize particles
void initParticles(fluid_particle **fluid_particle_pointers, fluid_particle *fluid_particles,
                   AABB* water, int start_x, int number_particles_x,    
                   edge *edges, int max_fluid_particles_local, float spacing, param* params)
{
    int i;
    fluid_particle *p;

    // Create fluid volume
    constructFluidVolume(fluid_particle_pointers, fluid_particles, water, start_x, number_particles_x, edges, spacing, params);

    // NULL out unused fluid pointers
    for(i=params->number_fluid_particles_local; i<max_fluid_particles_local; i++)
        fluid_particle_pointers[i] = NULL;

    // Initialize particle values
    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particle_pointers[i]->a_x = 0.0f;
        fluid_particle_pointers[i]->a_y = 0.0f;
        fluid_particle_pointers[i]->v_x = 0.0f;
        fluid_particle_pointers[i]->v_y = 0.0f;
    }
}
