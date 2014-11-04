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

#ifndef fluid_communication_h
#define fluid_communication_h

#include "structs.h"
#include "mpi.h"

// MPI globals
//MPI_Datatype Particletype;
MPI_Datatype TunableParamtype;
MPI_Datatype LeftEdgetype;
MPI_Datatype RightEdgetype;
MPI_Comm MPI_COMM_COMPUTE;
MPI_Group group_world;
MPI_Group group_compute;
MPI_Group group_render;

void create_MPI_types();
void create_communicators();
void free_MPI_types();
void halo_exchange(fluid_sim_t *fluid_sim);
void pack_halo_components(float *packed_recv_left, float *packed_recv_right, fluid_sim_t *fluid_sim);
void unpack_halo_components(float *packed_recv_left, float *packed_recv_right, fluid_sim_t *fluid_sim);
void pack_oob_components(float *left_send, float *right_send, fluid_sim_t *fluid_sim);
void unpack_oob_components(float *packed_recv, int num_recv, fluid_sim_t *fluid_sim);
void transfer_OOB_particles(fluid_sim_t *fluid_sim);
void update_halo_lambdas(fluid_sim_t *fluid_sim);
void update_halo_positions(fluid_sim_t *fluid_sim);

#endif
