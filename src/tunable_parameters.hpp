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

#ifndef tunable_parameters_h
#define tunable_parameters_h

#include <vector>
#include "renderer.hpp"

// enum of displayed parameter values
enum selected_param_t {
    MIN = 0,
    GRAVITY = MIN,
    SMOOTH,
    DENSITY,
    K,
    DQ,
    VISCOSITY,
    MAX = VISCOSITY
};

class TunableParameters
{
    public:
        TunableParameters(int num_compute_procs): num_compute_procs(num_compute_procs),  
                                                  num_compute_procs_active(num_compute_procs),
                                                  selected_parameter((selected_param_t)0){};
        ~TunableParameters();

        // Not quite sure how to deal with this...public for now
        selected_param_t selected_parameter;
        float rest_density;
        float smoothing_radius;
        float g;
        float k;
        float dq;
        float c;
        float time_step;
        float node_start_x;
        float node_end_x;
        float mover_center_x;
        float mover_center_y;
        float mover_center_z;
        float mover_radius;
        bool kill_sim;


        void move_parameter_up();
        void move_parameter_down();
        void increase_parameter();
        void decrease_parameter();
        void increase_gravity();
        void decrease_gravity();
        void increase_density();
        void decrease_density();
        void increase_viscosity();
        void decrease_viscosity();
        void increase_k();
        void decrease_k();
        void increase_dq();
        void decrease_dq();
        void increase_smoothing_radius();
        void decrease_smoothing_radius();
        void increase_mover_radius();
        void decrease_mover_radius();
        void add_partition();
        void remove_partition();
        void toggle_pause();
        void set_mover_gl_center(const float ogl_x, const float ogl_y, const float ogl_z);
        void reset_mover_size();
        void check_partition_left(int *particle_counts, int total_particles);

    private:
        int num_compute_procs;
        int num_compute_procs_active;
        std::vector<float> proc_starts;
        std::vector<float> proc_ends;
};

#endif
