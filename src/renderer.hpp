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

#ifndef fluid_renderer_h
#define fluid_renderer_h

#include "fluid.h"
#include "gl.hpp"
#include "camera.hpp"
#include "tunable_parameters.hpp"

class Renderer
{
    public:
        Renderer(int num_compute_procs): pause(false), 
                                         view_controls(false), 
                                         num_compute_procs(num_compute_procs),
                                         gl(), 
                                         camera(gl),
                                         tunable_parameters(num_compute_procs){};
        void start_rendering();
        void opengl_to_sim(float x, float y, float z, float *sim_x, float *sim_y, float *sim_z);
        void sim_to_opengl(float x, float y, float z, float *gl_x, float *gl_y, float *gl_z);
        void update_node_params();
        void set_activity_time();
        bool input_is_active();
        void update_inactive_state();
        void enable_view_controls();
        void disable_view_controls();
        void zoom_in_view();
        void zoom_out_view();
        void set_view_angle(const float x_pos, const float y_pos);
        void move_in_view();
        void move_out_view();
        const int screen_width() const { return this->gl.get_screen_width(); };
        const int screen_height() const { return this->gl.get_screen_height(); }; 
        void param_struct_to_class();
        void param_class_to_struct();
    private:
        GL gl;
        Camera camera;
        TunableParameters tunable_parameters;

        float sim_width;
        float sim_height;
        float sim_depth;

        int num_compute_procs;

        // Struct vector used to deal with tunable param class to mpi sendable type
        std::vector<tunable_parameters_t> tunable_param_structs;
        bool view_controls; // When shift is held mouse controls view
        bool pause;
        double last_activity_time; // Used to determine if simulation is being used or not
};

#endif
