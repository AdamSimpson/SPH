#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "fluid.h"
#include "hash.h"
#include "fileio.h"

int main(int argc, char *argv[])
{
    param params;
    AABB water_volume;
    AABB boundary;
    
    // Boundary box
    boundary.min_x = 0.0;
    boundary.max_x = 0.5;
    boundary.min_y = 0.0;
    boundary.max_y = 0.5;
    boundary.min_z = 0.0;
    boundary.max_z = 0.5;
    
    // water volume
    water_volume.min_x = 0.05;
    water_volume.max_x = 0.4;
    water_volume.min_y = 0.05;
    water_volume.max_y = 0.4;
    water_volume.min_z = 0.05;
    water_volume.max_z = 0.3;
    
    // Simulation parameters
    params.number_fluid_particles = 1000;
    params.rest_density = 1000.0;
    params.g = 9.8;
    params.alpha = 0.01;
    params.surface_tension =  0.01;
    params.number_steps = 100;
    params.time_step = 0.01;    

    // Mass of each particle
    double volume = (water_volume.max_x - water_volume.min_x) * (water_volume.max_y - water_volume.min_y) * (water_volume.max_z - water_volume.min_z);
    params.mass_particle = params.rest_density * (volume/params.number_fluid_particles);
    
    // Cube calculated spacing
    params.spacing_particle = pow(volume/params.number_fluid_particles,1.0/3.0);
    
    // Smoothing radius, h
    params.smoothing_radius = params.spacing_particle;
    
    // Total number of particles
    params.number_particles = params.number_fluid_particles;
    
    // Number of steps before frame needs to be written for 30 fps
    int steps_per_frame = (int)(1.0/(params.time_step*30.0));
    
    // Calculate speed of sound for simulation
    double max_height = water_volume.max_y;
    double max_velocity = sqrt(2.0*params.g*max_height);
    params.speed_sound = max_velocity/sqrt(0.01);

    // Minimum stepsize from Courant-Friedrichs-Lewy condition
    double recomend_step = 0.4 * params.smoothing_radius / (params.speed_sound * (1+ 0.6*params.alpha));
    printf("Using time step: %f, Minimum recomended %f\n",params.time_step, recomend_step);    

    // Allocate fluid particles array
    fluid_particle *fluid_particles = (fluid_particle*) malloc(params.number_fluid_particles * sizeof(fluid_particle));
    // Allocate neighbors array
    neighbor *neighbors = (neighbor*) malloc(params.number_fluid_particles * sizeof(neighbor));

    // Allocate new hash with all values zeroed
    params.length_fluid_hash = params.number_fluid_particles;
    uint2 *fluid_hash = calloc(params.length_fluid_hash, sizeof(uint2));

    // +1 added because range begins at 0
    params.grid_size_x = ceil((boundary.max_x - boundary.min_x) / params.smoothing_radius) + 1;
    params.grid_size_y = ceil((boundary.max_y - boundary.min_y) / params.smoothing_radius) + 1;
    params.grid_size_z = ceil((boundary.max_z - boundary.min_z) / params.smoothing_radius) + 1;
    uint64_t grid_size = params.grid_size_x * params.grid_size_y * params.grid_size_z;
    printf("maxx %f, y:%f ,z: %f\n", (boundary.max_x - boundary.min_x),(boundary.max_y - boundary.min_y) ,(boundary.max_z - boundary.min_z) );
    printf("grid size %llu x:%d,y:%d,z:%d \n", grid_size,params.grid_size_x,params.grid_size_y,params.grid_size_z);
    // Allocate hash start/end arrays
    uint2 *fluid_hash_positions = malloc(grid_size * sizeof(uint2));

    // Initialize particles
    initParticles(fluid_particles, neighbors, fluid_hash, fluid_hash_positions, &water_volume, &params);

    // Print some parameters
    printf("Volume: %f, p_volume: %f, Mass: %f, h: %f, f_parts: %d \n", volume, volume/params.number_particles, params.mass_particle, params.smoothing_radius, params.number_fluid_particles);
    
    // Main loop
    int n;
    int fileNum=0;
    for(n=0; n<params.number_steps; n++) {

        hash_fluid(fluid_particles,  neighbors, fluid_hash, fluid_hash_positions,  &params);
        
        updatePressures(fluid_particles, neighbors, &params);
        
        updateAccelerations(fluid_particles, neighbors, &params);
        
        updatePositions(fluid_particles, &boundary, &params);
        
        if (n % steps_per_frame == 0)
	    printf("Step: %d\n", n);
//            writeFile(fluid_particles, fileNum++, &params);
    }
    
    // Release memory
    free(fluid_particles);
    free(neighbors);
    free(fluid_hash);
    
    return 0;
}

////////////////////////////////////////////////////////////////////////////
// Smooting kernels
// http://www.fas.org/sgp/othergov/doe/lanl/lib-www/la-pubs/00538238.pdf
////////////////////////////////////////////////////////////////////////////

// B spline smoothing kernel
double W(fluid_particle* p, fluid_particle *q, double h)
{
    double r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y) + (p->z-q->z)*(p->z-q->z));
    double C = 1.0/(M_PI*h*h*h);
    double u = r/h;
    double val = 0.0;
    if(u >= 2.0)
        return val;
    else if(u < 1.0 )
        val = 1.0 - (3.0/2.0)*u*u + (3.0/4.0)*u*u*u;
    else if(u >= 1.0 && u < 2.0)
        val = (1.0/4.0) * pow(2.0-u,3.0);
    
    val *= C;
    return val;
}

// Gradient of B spline kernel
double del_W(fluid_particle* p, fluid_particle *q, double h)
{
    double r = sqrt((p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y) + (p->z-q->z)*(p->z-q->z));
    double C = 1.0/(M_PI * h*h*h);
    double u = r/h;
    double val = 0.0;
    if(u >= 2.0)
        return val;
    else if(u < 1.0 )
        val = -1.0/(h*h) * (3.0 - 9.0/4.0*u);
    else if(u >= 1.0 && u < 2.0)
        val = -3.0/(4.0*h*r) * pow(2.0-u,2.0);
    
    val *= C;
    return val;
}

// Project particle to AABB surface and reflect velocity
void reflectParticle(fluid_particle *p, param* params, double pen_depth, double *norm)
{
    double dt = params->time_step;
    double restitution = 0.0; // No slip condition, applicable to viscious fluids
    
    // project particle back onto surface
    p->x = p->x + pen_depth * norm[0];
    p->y = p->y + pen_depth * norm[1];
    p->z = p->z + pen_depth * norm[2];
    
    double v_mag = sqrt(p->v_x*p->v_x + p->v_y*p->v_y + p->v_z*p->v_z);
    double vDotn = p->v_x*norm[0] + p->v_y*norm[1] + p->v_z*norm[2];
    
    // Calculate new reflected velocity
    p->v_x = p->v_x - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[0]);
    p->v_y = p->v_y - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[1]);
    p->v_z = p->v_z - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[2]);
    
    p-> v_half_x = p-> v_half_x - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[0]);
    p-> v_half_y = p-> v_half_y - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[1]);
    p-> v_half_z = p-> v_half_z - ((1.0 + (restitution*pen_depth)/(dt*v_mag)) * vDotn * norm[2]);
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
     
        // Surface normal
        // Use local coordinates otherwise small round off errors from world->local frame throw off sgn calculation
        double n_x = (double)sgn(local_cp_x - local_x);
        double n_y = (double)sgn(local_cp_y - local_y);
        double n_z = (double)sgn(local_cp_z - local_z);
        double norm = sqrt(n_x*n_x + n_y*n_y + n_z*n_z);
        n_x = n_x/norm;
        n_y = n_y/norm;
        n_z = n_z/norm;

        double normal[3] = {n_x,n_y,n_z};
        reflectParticle(p, params, depth, normal);
    }
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
// http://www.control.auc.dk/~hempel/projects/pdf/sph1.pdf
// This makes initialization much more stable
////////////////////////////////////////////////////////////////////////////

double computeDensity(fluid_particle *p, fluid_particle *q, param *params)
{

    double v_x = (p->v_x - q->v_x);
    double v_y = (p->v_y - q->v_y);
    double v_z = (p->v_z - q->v_z);
    
    double density = params->mass_particle * del_W(p,q,params->smoothing_radius);
    double density_x = density * v_x * (p->x - q->x);
    double density_y = density * v_y * (p->y - q->y);
    double density_z = density * v_z * (p->z - q->z);
    
    density = (density_x + density_y + density_z)*params->time_step;


//    double density = params->mass_particle * W(p,q,params->smoothing_radius);   
    return density;
}

double computePressure(fluid_particle *p, param *params)
{
/*
    double gam = 7.0;
    double B = params->rest_density * params->speed_sound*params->speed_sound / gam;
    double pressure =  B * (pow((p->density/params->rest_density),gam) - 1.0);
        
    return pressure;
*/

    double pressure = 3* (p->density - params->rest_density);
    return pressure;
}

void updatePressures(fluid_particle *fluid_particles, neighbor *neighbors, param *params)
{
    int i, j;
    fluid_particle *p, *q;
    neighbor* n;
    
    for(i=0; i<params->number_fluid_particles; i++) {
        p = &fluid_particles[i];
        n = &neighbors[i];
//        p->density = 0;       
        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
            p->density += computeDensity(p,q,params);
        }
        p->pressure = computePressure(p,params);
    }
}

// Compute force(accleration) on fluid particle p by fluid particle q
void computeAcceleration(fluid_particle *p, fluid_particle *q, param *params)
{
    double accel;
    double h = params->smoothing_radius;
    
        // Pressure force
        accel = (p->pressure/(p->density*p->density) + q->pressure/(q->density*q->density)) * params->mass_particle * del_W(p,q,h);
        p->a_x -= accel * (p->x - q->x);
        p->a_y -= accel * (p->y - q->y);
        p->a_z -= accel * (p->z - q->z);
        
        // Viscosity force
        double VdotR = (p->v_x-q->v_x)*(p->x-q->x) + (p->v_y-q->v_y)*(p->y-q->y) + (p->v_z-q->v_z)*(p->z-q->z);
        if(VdotR < 0.0)
        {
            double nu = 2.0 * params->alpha * h * params->speed_sound / (p->density + q->density);
            double r2 = (p->x-q->x)*(p->x-q->x) + (p->y-q->y)*(p->y-q->y) + (p->z-q->z)*(p->z-q->z);
            double eps = params->spacing_particle/10.0; 
            double stress = nu * VdotR / (r2 + eps*h*h);
            accel = params->mass_particle * stress * del_W(p, q, h);
            p->a_x += accel * (p->x - q->x);
            p->a_y += accel * (p->y - q->y);
            p->a_z += accel * (p->z - q->z);
        }

        //Surface tension
        // BT 07 http://cg.informatik.uni-freiburg.de/publications/2011_GRAPP_airBubbles.pdf
        accel = params->surface_tension * W(p,q,h);
        p->a_x += accel * (p->x - q->x);
        p->a_y += accel * (p->y - q->y);
        p->a_z += accel * (p->z - q->z);

}

// Update particle acclerations
void updateAccelerations(fluid_particle *fluid_particles, neighbor *neighbors, param *params)
{
    int i,j;
    fluid_particle *p, *q;
    neighbor *n;
    
    for(i=0; i<params->number_fluid_particles; i++) {
        p = &fluid_particles[i];
        n = &neighbors[i];
        
        p->a_x = 0.0;
        p->a_y = 0.0;
        p->a_z = -params->g;
        
        // Acceleration on p due to neighbor fluid particles
        for(j=0; j<n->number_fluid_neighbors; j++) {
                q = n->fluid_neighbors[j];
                if (p!=q)
                    computeAcceleration(p,q,params);
        }
        
    }
}

// Leap Frog integration with v(t+1) estimated
void updateParticle(fluid_particle *p, param *params)
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
void updatePositions(fluid_particle *fluid_particles, AABB *boundary, param *params)
{
    int i;
    fluid_particle *p;
    
    for(i=0; i<params->number_fluid_particles; i++) {
        p = &fluid_particles[i];
        updateParticle(p, params);
	boundaryConditions(p, boundary, params);
    }
}

// Seed simulation with Euler step v(t-dt/2) needed by leap frog integrator
void eulerStart(fluid_particle* fluid_particles, neighbor* neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, param *params)
{
    //Generate neighbor hash
    hash_fluid(fluid_particles, neighbors, fluid_hash, fluid_hash_positions, params);
    
    updatePressures(fluid_particles, neighbors, params);
    
    updateAccelerations(fluid_particles, neighbors, params);
    
    // Set V (t0 - dt/2)
    int i;
    double dt_half;
    fluid_particle *p;
    
    for(i=0; i<params->number_fluid_particles; i++)
    {
        p = &fluid_particles[i];
        dt_half = params->time_step/2.0;
        p->v_half_x = p->v_x - p->a_x * dt_half;
        p->v_half_y = p->v_y - p->a_y * dt_half;
        p->v_half_z = p->v_z - p->a_z * dt_half;
    }
}

// Initialize particles
void initParticles(fluid_particle* fluid_particles, neighbor *neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, AABB* water, param* params)
{
    double spacing = params->spacing_particle;

    // Initialize particle values
    int i;
    for(i=0; i<params->number_fluid_particles; i++) {
        fluid_particles[i].a_x = 0.0;
        fluid_particles[i].a_y = 0.0;
        fluid_particles[i].a_z = 0.0;
        fluid_particles[i].v_x = 0.0;
        fluid_particles[i].v_y = 0.0;
        fluid_particles[i].v_z = 0.0;
        fluid_particles[i].density = params->rest_density;
    }

    // Place particles inside bounding volume
    double x,y,z;
    i = 0;
    for(z=water->min_z; z<=water->max_z; z+=spacing) {
        for(y=water->min_y; y<=water->max_y; y+=spacing) {
            for(x=water->min_x; x<=water->max_x; x+=spacing) {
                if(i < params->number_fluid_particles) {
                    fluid_particles[i].x = x;
                    fluid_particles[i].y = y;
                    fluid_particles[i].z = z;
                    i++;
                }
            }
        }
    }
    params->number_fluid_particles = i;
    
    // Set initial velocity half step
    eulerStart(fluid_particles, neighbors, fluid_hash, fluid_hash_positions, params);
}

////////////////////////////////////////////////
// Utility Functions
////////////////////////////////////////////////
double min(double a, double b){
    double min = a;
    min = b < min ? b : min;
    return min;
}

double max(double a, double b){
    double max = a;
    max = b > max ? b : max;
    return max;
}

int sgn(double x) {
    int val = 0;
    if (x < 0.0)
	val = -1;
    else if (x > 0.0)
	val = 1;

    return val;
}
