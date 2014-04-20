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

#ifndef controls_h
#define controls_h

#include "renderer.h"

void move_parameter_up(render_t *render_state);
void move_parameter_down(render_t *render_state);
void increase_parameter(render_t *render_state);
void decrease_parameter(render_t *render_state);
void increase_gravity(render_t *render_state);
void decrease_gravity(render_t *render_state);
void increase_density(render_t *render_state);
void decrease_density(render_t *render_state);
void increase_viscosity(render_t *render_state);
void decrease_viscosity(render_t *render_state);
void increase_pressure(render_t *render_state);
void decrease_pressure(render_t *render_state);
void increase_elasticity(render_t *render_state);
void decrease_elasticity(render_t *render_state);
void increase_mover_width(render_t *render_state);
void decrease_mover_width(render_t *render_state);
void increase_mover_height(render_t *render_state);
void decrease_mover_height(render_t *render_state);
void set_fluid_x(render_t *render_state);
void set_fluid_y(render_t *render_state);
void set_fluid_a(render_t *render_state);
void set_fluid_b(render_t *render_state);
void remove_partition(render_t *render_state);
void add_partition(render_t *render_state);
void toggle_dividers(render_t *state);
void toggle_pause(render_t *state);
void set_mover_gl_center(render_t *render_state, float ogl_x, float ogl_y);


#endif