#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fluid.h"


// Compiler must support vardiac macros
#ifdef DEBUG
# define DEBUG_PRINT(str, args...) printf("DEBUG: %s:%d:%s(): " str, \
 __FILE__, __LINE__, __func__, ##args)
#else
# define DEBUG_PRINT(str, args...) do {} while (0)
#endif

////////////////////////////////////////////////////////////////////////////
// Smoothing Kernels
// Don't use pow as it's POWerfully slow
///////////////////////////////////////////////////////////////////////////

// (h^2 - r^2)^3 normalized in 3D (poly6)
double W(double r, double h)
{
    if(r > h)
        return 0.0;

    double C = 315.0/(64.0*M_PI* h*h*h*h*h*h*h*h*h);
    double W = C*(h*h-r*r)*(h*h-r*r)*(h*h-r*r);
    return W;
}

// Gradient (h-r)^3 normalized in 3D (Spikey) magnitude
// Need to multiply by r/|r|
double del_W(double r, double h)
{
    if(r > h)
        return 0.0;

    double C = -45.0/(M_PI * h*h*h*h*h*h);
    double del_W = C*(h-r)*(h-r);
    return del_W;
}

////////////////////////////////////////////////////////////////////////////
// Particle attribute computations
////////////////////////////////////////////////////////////////////////////

void vorticity_confinement(fluid_particle_t *fluid_particles, neighbor_t* neighbors, param_t *params)
{
    int i,j;
    fluid_particle_t *p, *q;
    neighbor_t *n;
    double dt = params->time_step;
    double x_diff, y_diff, z_diff, vx_diff, vy_diff, vz_diff,
          r_mag, dw, dw_x, dw_y, dw_z, part_vort_x, part_vort_y, part_vort_z,
          vort_x, vort_y, vort_z, eta_x, eta_y, eta_z, eta_mag, N_x, N_y, N_z;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        p = &fluid_particles[i];
        n = &neighbors[i];
        vort_x = 0.0;
        vort_y = 0.0;
        vort_z = 0.0;
        eta_x  = 0.0;
        eta_y  = 0.0;
        eta_z  = 0.0;

        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q = &fluid_particles[n->neighbor_indices[j]];

            x_diff = p->x_star - q->x_star;
            y_diff = p->y_star - q->y_star;
            z_diff = p->z_star - q->z_star;

            vx_diff = q->v_x - p->v_x;
            vy_diff = q->v_y - p->v_y;
            vz_diff = q->v_z - p->v_z;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);

            dw = del_W(r_mag, params->smoothing_radius);
            if(r_mag < 0.0001) {
              printf("p->x_star: %f\n", p->x_star);
              r_mag = 0.0001;
            }
            dw_x = dw*x_diff/r_mag;
            dw_y = dw*y_diff/r_mag;
            dw_z = dw*z_diff/r_mag;

            part_vort_x =  vy_diff*dw_z - vz_diff*dw_y;
            part_vort_y =  vz_diff*dw_x - vx_diff*dw_z;
            part_vort_z =  vx_diff*dw_y - vy_diff*dw_x;

            vort_x += part_vort_x;
            vort_y += part_vort_y;
            vort_z += part_vort_z;

            if(x_diff<0.0000001 || y_diff<0.0000001 || z_diff<0.0000001)
                continue;

            eta_x += fabs(part_vort_x)/x_diff;
            eta_y += fabs(part_vort_y)/y_diff;
            eta_z += fabs(part_vort_z)/z_diff;
        }
        eta_mag = sqrt(eta_x*eta_x + eta_y*eta_y + eta_z*eta_z);
        if(eta_mag<0.0000001)
            continue;

        N_x = eta_x / eta_mag;
        N_y = eta_y / eta_mag;
        N_z = eta_z / eta_mag;
        p->v_x += dt * ( N_y*vort_z - N_z*vort_y);
        p->v_y += dt * ( N_z*vort_x - N_x*vort_z);
        p->v_z += dt * ( N_x*vort_y - N_y*vort_x);
    }
}

void XSPH_viscosity(fluid_particle_t *fluid_particles, neighbor_t* neighbors, param_t *params)
{
    int i,j;
    unsigned int q_index;
    neighbor_t *n;
    double c = params->c;
    double h = params->smoothing_radius;

    double x_diff, y_diff, z_diff, vx_diff, vy_diff, vz_diff, r_mag, w;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        double partial_sum_x = 0.0;
        double partial_sum_y = 0.0;
        double partial_sum_z = 0.0;
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->neighbor_indices[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            vx_diff = fluid_particles[q_index].v_x - fluid_particles[i].v_x;
            vy_diff = fluid_particles[q_index].v_y - fluid_particles[i].v_y;
            vz_diff = fluid_particles[q_index].v_z - fluid_particles[i].v_z;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff * z_diff*z_diff);

            w = W(r_mag, h);

            // http://mmacklin.com/pbf_sig_preprint.pdf is missing 1/sigma contribution
            // see: http://www.cs.ubc.ca/~rbridson/docs/schechter-siggraph2012-ghostsph.pdf
            partial_sum_x += vx_diff * w / fluid_particles[q_index].density;
            partial_sum_y += vy_diff * w / fluid_particles[q_index].density;
            partial_sum_z += vz_diff * w / fluid_particles[q_index].density;
        }

        partial_sum_x *= c;
        partial_sum_y *= c;
        partial_sum_z *= c;

        fluid_particles[i].v_x += partial_sum_x;
        fluid_particles[i].v_y += partial_sum_y;
        fluid_particles[i].v_z += partial_sum_z;
    }
}

void compute_densities(fluid_particle_t *fluid_particles, neighbor_t *neighbors, param_t *params)
{
    int i,j;
    unsigned int q_index;
    neighbor_t *n;
    double h = params->smoothing_radius;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        double x_diff, y_diff, z_diff, r_mag, density;

        // Own contribution to density
        density = W(0.0, h);

        // Neighbor contribution
        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->neighbor_indices[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            density += W(r_mag, h);
        }

        // Update particle density
        fluid_particles[i].density = density;
    }
}

void apply_gravity(fluid_particle_t *fluid_particles, param_t *params)
{
    int i;
    double dt = params->time_step;
    double g = -params->g;

    for(i=0; i<(params->number_fluid_particles_local); i++) {
        fluid_particles[i].v_y += g*dt;
    }
}

void update_dp_positions(fluid_particle_t *fluid_particles, AABB_t *boundary_global, param_t *params)
{
    int i;

    for(i=0; i<(params->number_fluid_particles_local); i++) {

        fluid_particles[i].x_star += fluid_particles[i].dp_x;
        fluid_particles[i].y_star += fluid_particles[i].dp_y;
        fluid_particles[i].z_star += fluid_particles[i].dp_z;

        // Enforce boundary conditions
        boundary_conditions(fluid_particles, i, boundary_global);

    }
}

void update_positions(fluid_particle_t *fluid_particles, param_t *params)
{
     int i;

     for(i=0; i<(params->number_fluid_particles_local); i++) {
        fluid_particles[i].x = fluid_particles[i].x_star;
        fluid_particles[i].y = fluid_particles[i].y_star;
        fluid_particles[i].z = fluid_particles[i].z_star;
    }
}

void calculate_lambda(fluid_particle_t *fluid_particles, neighbor_t *neighbors, param_t *params)
{
    int i,j;
    unsigned int q_index;
    neighbor_t *n;

    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        double Ci = fluid_particles[i].density/params->rest_density - 1.0;

        double sum_C, x_diff, y_diff, z_diff, r_mag,
              grad, grad_x, grad_y, grad_z,
              sum_grad_x, sum_grad_y, sum_grad_z;

        sum_C = 0.0;
        grad_x = 0.0;
        grad_y = 0.0;
        grad_z = 0.0;
        sum_grad_x = 0.0;
        sum_grad_y = 0.0;
        sum_grad_z = 0.0;

        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->neighbor_indices[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);

            grad = del_W(r_mag, params->smoothing_radius);
            if(r_mag < 0.0001) {
              r_mag = 0.0001;
              printf("pstar: %f, %f, %f grad: %f\n", fluid_particles[i].x_star, fluid_particles[i].y_star,fluid_particles[i].z_star, grad);
            }
            grad_x = grad*x_diff/r_mag;
            grad_y = grad*y_diff/r_mag;
            grad_z = grad*z_diff/r_mag;
            sum_grad_x += grad_x;
            sum_grad_y += grad_y;
            sum_grad_z += grad_z;
            // Add k = j contribution
            sum_C += (grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
        }

        // k = i contribution
        sum_C += sum_grad_x*sum_grad_x
              + sum_grad_y*sum_grad_y
              + sum_grad_z*sum_grad_z;

        sum_C *= (1.0/(params->rest_density*params->rest_density));

        double epsilon = 1.0;
        fluid_particles[i].lambda = -Ci/(sum_C + epsilon);
    }
}

void update_dp(fluid_particle_t *fluid_particles, neighbor_t *neighbors, param_t *params)
{
    unsigned int q_index;
    neighbor_t *n;
    double x_diff, y_diff, z_diff, dp, r_mag;
    double h = params->smoothing_radius;
    double k = params->k;
    double dq = params->dq;
    double Wdq = W(dq, h);

    int i,j;
    for(i=0; i<params->number_fluid_particles_local; i++)
    {
        n = &neighbors[i];

        double dp_x = 0.0;
        double dp_y = 0.0;
        double dp_z = 0.0;
        double s_corr;
        double WdWdq; // W()/Wdq()

        for(j=0; j<n->number_fluid_neighbors; j++)
        {
            q_index = n->neighbor_indices[j];
            x_diff = fluid_particles[i].x_star - fluid_particles[q_index].x_star;
            y_diff = fluid_particles[i].y_star - fluid_particles[q_index].y_star;
            z_diff = fluid_particles[i].z_star - fluid_particles[q_index].z_star;

            r_mag = sqrt(x_diff*x_diff + y_diff*y_diff + z_diff*z_diff);
            WdWdq = W(r_mag, h)/Wdq;
            s_corr = -k*WdWdq*WdWdq*WdWdq*WdWdq;
            dp = (fluid_particles[i].lambda + fluid_particles[q_index].lambda + s_corr)*del_W(r_mag, h);

            if(r_mag < 0.0001) {
              printf("pstar: %f, %f, %f grad: %f\n", fluid_particles[i].x_star, fluid_particles[i].y_star,fluid_particles[i].z_star, WdWdq);
              r_mag = 0.0001;
            }
            dp_x += dp*x_diff/r_mag;
            dp_y += dp*y_diff/r_mag;
            dp_z += dp*z_diff/r_mag;
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
        if (fluid_particles[i].x_star < params->node_start_x)
            out_of_bounds->oob_indices_left[out_of_bounds->number_oob_particles_left++] = i;
        else if (fluid_particles[i].x_star > params->node_end_x)
            out_of_bounds->oob_indices_right[out_of_bounds->number_oob_particles_right++] = i;
    }

   // Transfer particles that have left the processor bounds
   transferOOBParticles(fluid_particles, out_of_bounds, params);
}

// Predict position
void predict_positions(fluid_particle_t *fluid_particles, AABB_t *boundary_global, param_t *params)
{
    int i;
    double dt = params->time_step;

    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particles[i].x_star = fluid_particles[i].x + (fluid_particles[i].v_x * dt);
        fluid_particles[i].y_star = fluid_particles[i].y + (fluid_particles[i].v_y * dt);
        fluid_particles[i].z_star = fluid_particles[i].z + (fluid_particles[i].v_z * dt);

        // Enforce boundary conditions before hash
        // Otherwise predicted position can blow up hash
        boundary_conditions(fluid_particles, i, boundary_global);
    }
}

void check_velocity(double *v_x, double *v_y, double *v_z)
{
    double v_max = 30.0;

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
    int i;
    double dt = params->time_step;
    double v_x, v_y, v_z;

    // Update local and halo particles, update halo so that XSPH visc. is correct
    for(i=0; i<params->number_fluid_particles_local + params->number_halo_particles_left + params->number_halo_particles_right; i++) {
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
void boundary_conditions(fluid_particle_t *fluid_particles, unsigned int i, AABB_t *boundary)
{

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(fluid_particles[i].x_star  < boundary->min_x)
        fluid_particles[i].x_star = boundary->min_x;
    else if(fluid_particles[i].x_star  > boundary->max_x)
        fluid_particles[i].x_star = boundary->max_x-0.00001;

    if(fluid_particles[i].y_star  <  boundary->min_y)
        fluid_particles[i].y_star = boundary->min_y;
    else if(fluid_particles[i].y_star  > boundary->max_y)
        fluid_particles[i].y_star = boundary->max_y-0.00001;

    if(fluid_particles[i].z_star  <  boundary->min_z)
        fluid_particles[i].z_star = boundary->min_z;
    else if(fluid_particles[i].z_star  > boundary->max_z)
        fluid_particles[i].z_star = boundary->max_z-0.00001;
}

// Initialize particles
void initParticles(fluid_particle_t *fluid_particles,
                neighbor_t *neighbors, bucket_t *hash, AABB_t* water, AABB_t* boundary_global,
                int start_x, int number_particles_x, edge_t *edges, param_t* params)
{
    int i;

    // Create fluid volume
    constructFluidVolume(fluid_particles, water, start_x, number_particles_x, edges, params);

    // Initialize particle values
    for(i=0; i<params->number_fluid_particles_local; i++) {
        fluid_particles[i].x_star = fluid_particles[i].x;
        fluid_particles[i].y_star = fluid_particles[i].y;
        fluid_particles[i].z_star = fluid_particles[i].z;
        fluid_particles[i].v_x = 0.0;
        fluid_particles[i].v_y = 0.0;
        fluid_particles[i].v_z = 0.0;
        fluid_particles[i].dp_x = 0.0;
        fluid_particles[i].dp_y = 0.0;
        fluid_particles[i].dp_z = 0.0;
        fluid_particles[i].lambda = 0.0;
        fluid_particles[i].density = params->rest_density;
    }
}
