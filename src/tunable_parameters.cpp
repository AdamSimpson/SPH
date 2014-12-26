/*
The MIT License (MIT)

Copyright (c) 2014 Adam Simpson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include "tunable_parameters.hpp"
#include "renderer.hpp"

// Move selected parameter up
void TunableParameters::move_parameter_up()
{
    if(this->selected_parameter == MIN)
        this->selected_parameter = MAX;
    else
        this->selected_parameter = (selected_param_t)((int)this->selected_parameter-1);
}

// Move selected parameter down
void TunableParameters::move_parameter_down() 
{
    if(this->selected_parameter == MAX)
        this->selected_parameter = MIN;
    else
        this->selected_parameter = (selected_param_t)((int)this->selected_parameter+1);
}

void TunableParameters::increase_parameter()
{
    switch(this->selected_parameter) {
        case GRAVITY:
            this->increase_gravity();
            break;
        case SMOOTH:
            this->increase_smoothing_radius();
            break;
        case DENSITY:
            this->increase_density();
            break;
        case K:
            this->increase_k();
            break;
        case DQ:
            this->increase_dq();
            break;
         case VISCOSITY:
            this->increase_viscosity();
            break;
    }
}

void TunableParameters::decrease_parameter()
{
    switch(this->selected_parameter) {
        case GRAVITY:
            this->decrease_gravity();
            break;
        case SMOOTH:
            this->decrease_smoothing_radius();
            break;
        case DENSITY:
            this->decrease_density();
            break;
        case K:
            this->decrease_k();
            break;
        case DQ:
            this->decrease_dq();
            break;
        case VISCOSITY:
            this->decrease_viscosity();
            break;
    }
}

// Increase gravity parameter
void TunableParameters::increase_gravity()
{
    static const float max_grav = -9.0f;
    if(this->g > max_grav)
        this->g -= 1.0f;
}

// Decreate gravity parameter
void TunableParameters::decrease_gravity()
{
    static const float min_grav = 9.0f;
    if(this->g < min_grav)
        this->g -= 1.0f;
}

// Increase density parameter
void TunableParameters::increase_density()
{
    static const float max_dens = 5.0f;
    if(this->rest_density < max_dens)
        this->rest_density += 0.01f;
}

// Decreate gravity parameter
void TunableParameters::decrease_density()
{
    static const float min_dens = -5.0f;
    if(this->rest_density > min_dens)
        this->rest_density -= 0.01f;
}

// Decreate smoothing_radius parameter
void TunableParameters::decrease_smoothing_radius()
{
    float min_radius = 0.0f;
    if(this->smoothing_radius > min_radius)
        this->smoothing_radius -= 0.1f;
}

// Increase smoothing_radius parameter
void TunableParameters::increase_smoothing_radius()
{
    float max_radius = 5.0f;
    if(this->smoothing_radius < max_radius)
        this->smoothing_radius += 0.1f;
}

// Decreate dq parameter
// Deal with multiples of the smoothing radius
void TunableParameters::decrease_dq()
{
    float min_dq = 0.0f;
    if(this->dq > min_dq)
        this->dq -= (0.05f*this->smoothing_radius);
}

// Increase dq parameter
// Deal with multiples of the smoothing radius
void TunableParameters::increase_dq()
{
    float max_dq = 1.0f*this->smoothing_radius;
    if(dq < max_dq)
        this->dq += (0.05f*this->smoothing_radius);
}


// Increase viscosity parameter
void TunableParameters::increase_viscosity()
{
    static const float max_viscosity = 100.0f;
    if(this->c < max_viscosity)
        this->c += 0.05f;
}

// Decreate viscosity parameter
void TunableParameters::decrease_viscosity()
{
    static const float min_viscosity = -100.0f;
    if(this->c > min_viscosity)
        this->c -= 0.05f;
}

// Increase k parameter
void TunableParameters::increase_k()
{
    static const float max_k = 5.0f;
    if(this->k < max_k)
        this->k += 0.05f;
}

// Decreate k parameter
void TunableParameters::decrease_k()
{
    static const float min_k = -5.0f;
    if(this->k > min_k)
        this->k -= 0.05f;
}

// Set center of mover, input is openGL coordinates
void TunableParameters::set_mover_gl_center(const float ogl_x, const float ogl_y, const float ogl_z)
{
    float sim_x, sim_y, sim_z;
//    opengl_to_sim(ogl_x, ogl_y, ogl_z, &sim_x, &sim_y, &sim_z);

    this->mover_center_x = sim_x;
    this->mover_center_y = sim_y;
    this->mover_center_z = sim_z;
}

// Increase mover radius
void TunableParameters::increase_mover_radius()
{
    // Maximum radius of mover
    static const float max_width = 4.0f;
    if(this->mover_radius < max_width)
        this->mover_radius += 0.2f;
}

// Decrease mover radius
void TunableParameters::decrease_mover_radius()
{
    // Minimum width of mover
    static const float min_width = 1.0f;
    if(this->mover_radius > min_width)
        this->mover_radius -= 0.2f;
}

// Reset the mover radius
void TunableParameters::reset_mover_radius()
{
    this->mover_radius = 2.0f;
}

// Set last partition to be outside of simulation bounds
// Effectively removing it from the simulation
void TunableParameters::remove_partition()
{
    if(this->num_compute_procs_active == 1) 
	return;

    int removed_rank = this->num_compute_procs_active-1;

    // Set new end position of last active proc to end of simulation
    this->proc_ends[removed_rank-1] = this->proc_ends[removed_rank];

    // Send start and end x out of sim bounds for removed proc
    float position = this->proc_ends[removed_rank] + 1.0; // +1.0 ensures it's out of the simulation bounds
    this->proc_starts[removed_rank] = position;
    this->proc_ends[removed_rank] = position;

    this->num_compute_procs_active -= 1;
}

// Add on partition to right side that has been removed
void TunableParameters::add_partition()
{

    if(this->num_compute_procs_active == this->num_compute_procs)
	return;

    // Length of currently last partiion
    float length = this->proc_ends[this->num_compute_procs_active-1] - this->proc_starts[this->num_compute_procs_active-1];

    // If the last partition is too small we can't split it into two
    if(length < 2.5*this->smoothing_radius)
	return;

    // Set end of added partition to current end location
    this->proc_ends[this->num_compute_procs_active] = this->proc_ends[this->num_compute_procs_active-1];
    
    // Divide the current last partition in half
    float new_x = proc_starts[this->num_compute_procs_active-1] + length*0.5;
    this->proc_ends[this->num_compute_procs_active-1] = new_x;
    this->proc_starts[this->num_compute_procs_active] = new_x;

    this->num_compute_procs_active += 1;
}

// Checks for a balanced number of particles on each compute node
// If unbalanced the partition between nodes will change 
// Check from right to left
void TunableParameters::check_partition_left(int *particle_counts, int total_particles)
{
    int rank, diff;
    float h, dx, length, length_left, length_right;

    // Particles per proc if evenly divided
    int even_particles = total_particles/this->num_compute_procs_active;
    int max_diff = even_particles/15.0f;

    // Fixed distance to move partition is 0.125*smoothing radius
    h = this->smoothing_radius;
    dx = h*0.125;

    for(rank=this->num_compute_procs_active; rank-- > 1; )
    {
        length =  this->proc_ends[rank] - this->proc_starts[rank];
        length_left =  this->proc_ends[rank-1] - this->proc_starts[rank-1];
        diff = particle_counts[rank] - even_particles;

        // current rank has too many particles
        if( diff > max_diff && length > 2*h) {
            this->proc_starts[rank] += dx;
            this->proc_ends[rank-1] = this->proc_starts[rank];
        }
        // current rank has too few particles
        else if (diff < -max_diff && length_left > 2*h) {
            this->proc_starts[rank] -= dx;
            this->proc_ends[rank-1] = this->proc_starts[rank];
        }
    }

    // Left most partition has a tendency to not balance correctly so we test it explicitly
    if(this->num_compute_procs_active > 1)
    {
        length =  this->proc_ends[0] - this->proc_starts[0];
        length_right =  this->proc_ends[1]- this->proc_starts[1];
        diff = particle_counts[0] - even_particles;

        // current rank has too many particles
        if( diff > max_diff && length > 2*h) {
            this->proc_ends[0] -= dx;
            this->proc_starts[1] = this->proc_ends[0];
        }
        else if (diff < -max_diff && length_right > 2*h) {
            this->proc_ends[0] += dx;
            this->proc_starts[1] = this->proc_ends[0];
        }
    }
}
