#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include "fluid.h"
#include "hash.h"
#include "geometry.h"
#include "fileio.h"

int main(int argc, char *argv[])
{
    param params;
    AABB water_volume;
    AABB boundary;
    
    // Boundary box
    boundary.min_x = 0.0;
    boundary.max_x = 1.1;
    boundary.min_y = 0.0;
    boundary.max_y = 1.1;
    boundary.min_z = 0.0;
    boundary.max_z = 1.1;
    
    // water volume
    water_volume.min_x = 0.1;
    water_volume.max_x = 0.5;
    water_volume.min_y = 0.1;
    water_volume.max_y = 0.5;
    water_volume.min_z = 0.1;
    water_volume.max_z = 0.8;
    
    // Simulation parameters
    params.number_fluid_particles = 1600;
    params.rest_density = 1000.0;
    params.g = 9.8;
    params.alpha = 0.02;
    params.surface_tension =  0.01;
    params.number_steps = 1000;
    params.time_step = 0.00015;    

    // Mass of each particle
    double volume = (water_volume.max_x - water_volume.min_x) * (water_volume.max_y - water_volume.min_y) * (water_volume.max_z - water_volume.min_z);
    params.mass_particle = params.rest_density * (volume/params.number_fluid_particles);
    
    // Cube calculated spacing
    params.spacing_particle = pow(volume/params.number_fluid_particles,1.0/3.0);
    
    // Smoothing radius, h
    params.smoothing_radius = params.spacing_particle;
    
    // Boundary particles
    int num_x = ceil((boundary.max_x - boundary.min_x)/params.spacing_particle);
    int num_y = ceil((boundary.max_y - boundary.min_y)/params.spacing_particle);
    int num_z = ceil((boundary.max_z - boundary.min_z)/params.spacing_particle);
    int num_boundary_particles = (2 * num_x * num_z) + (2 * num_y * num_z) + (2* num_y * num_z);
    params.number_boundary_particles = num_boundary_particles;
    
    // Total number of particles
    params.number_particles = params.number_boundary_particles + params.number_fluid_particles;
    
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
    // Allocate boundary particles array
    boundary_particle *boundary_particles = (boundary_particle*) malloc(params.number_boundary_particles * sizeof(boundary_particle));
    // Allocate neighbors array
    neighbor *neighbors = (neighbor*) malloc(params.number_fluid_particles * sizeof(neighbor));

    // Allocate new hash with all values zeroed
    params.length_fluid_hash = params.number_fluid_particles;
    uint2 *fluid_hash = calloc(params.length_fluid_hash, sizeof(uint2));
    params.length_boundary_hash = params.number_boundary_particles;
    uint2 *boundary_hash = calloc(params.length_boundary_hash, sizeof(uint2));

    // Construct bounding box
    constructBoundaryBox(boundary_particles, &boundary, &params);

    // +1 added because range begins at 0
    params.grid_size_x = ceil((boundary.max_x - boundary.min_x) / params.smoothing_radius) + 1;
    params.grid_size_y = ceil((boundary.max_y - boundary.min_y) / params.smoothing_radius) + 1;
    params.grid_size_z = ceil((boundary.max_z - boundary.min_z) / params.smoothing_radius) + 1;
    uint64_t grid_size = params.grid_size_x * params.grid_size_y * params.grid_size_z;
    printf("maxx %f, y:%f ,z: %f\n", (boundary.max_x - boundary.min_x),(boundary.max_y - boundary.min_y) ,(boundary.max_z - boundary.min_z) );
    printf("grid size %llu x:%d,y:%d,z:%d \n", grid_size,params.grid_size_x,params.grid_size_y,params.grid_size_z);
    // Allocate hash start/end arrays
    uint2 *fluid_hash_positions = malloc(grid_size * sizeof(uint2));
    uint2 *boundary_hash_positions = malloc(grid_size * sizeof(uint2));

    // Initialize particles
    initParticles(fluid_particles, boundary_particles, neighbors, fluid_hash, fluid_hash_positions, boundary_hash, boundary_hash_positions, &water_volume, &boundary, &params);

    // Write boundary particles to file
    writeBoundaryFile(boundary_particles, &params);
    
    // Print some parameters
    printf("Volume: %f, p_volume: %f, Mass: %f, h: %f, b_parts: %d, f_parts: %d, t_parts %d\n", volume, volume/params.number_particles, params.mass_particle, params.smoothing_radius, params.number_boundary_particles, params.number_fluid_particles,params.number_particles);
    
    // Main loop
    int n;
    int fileNum=0;
    for(n=0; n<params.number_steps; n++) {

        hash_fluid(fluid_particles, boundary_particles,  neighbors, fluid_hash, fluid_hash_positions, boundary_hash, boundary_hash_positions,  &params);
        
        updatePressures(fluid_particles, neighbors, &params);
        
        updateAccelerations(fluid_particles, neighbors, &params);
        
        updatePositions(fluid_particles, &params);
        
        if (n % steps_per_frame == 0)
            writeFile(fluid_particles, fileNum++, &params);
    }
    
    // Release memory
    free(fluid_particles);
    free(boundary_particles);
    free(neighbors);
    free(fluid_hash);
    free(boundary_hash);
    
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

////////////////////////////////////////////////////////////////////////////
// Boundary particle force
// http://iopscience.iop.org/0034-4885/68/8/R01/pdf/0034-4885_68_8_R01.pdf
////////////////////////////////////////////////////////////////////////////

double boundaryGamma(fluid_particle *p, boundary_particle *k, param *params)
{
    // Radial distance between p,q
    double r = sqrt((p->x-k->x)*(p->x-k->x) + (p->y-k->y)*(p->y-k->y) + (p->z-k->z)*(p->z-k->z));
    // Distance to p normal to surface particle
    double y = sqrt((p->x-k->x)*(p->x-k->x)*(k->n_x*k->n_x) + (p->y-k->y)*(p->y-k->y)*(k->n_y*k->n_y) + (p->z-k->z)*(p->z-k->z)*(k->n_z*k->n_z));
    // Tangential distance
    double x = r-y;
    
    double u = y/params->smoothing_radius;
    double xi = (1-x/params->spacing_particle)?x<params->spacing_particle:0.0;
    double C = xi*2.0*0.02 * params->speed_sound * params->speed_sound / y;
    double val = 0.0;
    
    if(u > 0.0 && u < 2.0/3.0)
        val = 2.0/3.0;
    else if(u < 1.0 && u > 2.0/3.0 )
        val = (2*u - 3.0/2.0*u*u);
    else if (u < 2.0 && u > 1.0)
        val = 0.5*(2.0-u)*(2.0-u);
    else
        val = 0.0;
    
    val *= C;
    
    return val;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
// http://www.control.auc.dk/~hempel/projects/pdf/sph1.pdf
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
    
    return density;
}

double computePressure(fluid_particle *p, param *params)
{
    double gam = 7.0;
    double B = params->rest_density * params->speed_sound*params->speed_sound / gam;
    double pressure =  B * (pow((p->density/params->rest_density),gam) - 1.0);
        
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
        
        for(j=0; j<n->number_fluid_neighbors; j++) {
            q = n->fluid_neighbors[j];
            p->density += computeDensity(p,q,params);
        }
        p->pressure = computePressure(p,params);
    }
}

void computeBoundaryAcceleration(fluid_particle *p, boundary_particle *k, param *params)
{
    double bGamma = boundaryGamma(p,k,params);
    p->a_x += bGamma * k->n_x;
    p->a_y += bGamma * k->n_y;
    p->a_z += bGamma * k->n_z;
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
    boundary_particle *k;
    
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
        
        // Acceleration on p due to neighbor boundary particles
        for (j=0; j<n->number_boundary_neighbors; j++) {
            k = n->boundary_neighbors[j];
            computeBoundaryAcceleration(p,k,params);
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
void updatePositions(fluid_particle *fluid_particles, param *params)
{
    int i;
    fluid_particle *p;
    
    for(i=0; i<params->number_fluid_particles; i++) {
        p = &fluid_particles[i];
        updateParticle(p, params);
    }
}

// Seed simulation with Euler step v(t-dt/2) needed by leap frog integrator
void eulerStart(fluid_particle* fluid_particles, boundary_particle* boundary_particles, neighbor* neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions,  uint2 *boundary_hash, uint2 *boundary_hash_positions, param *params)
{
    //Generate neighbor hash
    hash_fluid(fluid_particles, boundary_particles, neighbors, fluid_hash, fluid_hash_positions, boundary_hash, boundary_hash_positions, params);
    
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
void initParticles(fluid_particle* fluid_particles, boundary_particle* boundary_particles,
                neighbor *neighbors, uint2 *fluid_hash, uint2 *fluid_hash_positions, uint2 *boundary_hash, uint2 *boundary_hash_positions, AABB* water, AABB* boundary, param* params)
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
    
    // Place boundary particles in hash
    hash_boundary(boundary_particles, boundary_hash, boundary_hash_positions, params);
    
    // Set initial velocity half step
    eulerStart(fluid_particles, boundary_particles, neighbors, fluid_hash, fluid_hash_positions, boundary_hash, boundary_hash_positions, params);
}
