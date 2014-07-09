#include "thrust/device_ptr.h"
#include "thrust/sort.h"
#include "fluid.h"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// Calculate the density contribution of p on q and q on p
__global__ void calculate_density(fluid_particle **fluid_particle_pointers, uint *start_indexes, uint *end_indexes, uint *particle_ids, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;

    int num_fluid, grid_x, grid_y, bucket_index;
    uint start_index, end_index;
    fluid_particle *p, *q;
    float p_x, p_y, ratio, QmP_x, QmP_y, OmR2, r, h_recip;
    h_recip = 1.0f/params->tunable_params.smoothing_radius;

    num_fluid = params->number_fluid_particles_local + params->number_halo_particles;

    if(i >= num_fluid)
        return;

    p = fluid_particle_pointers[i];
    p_x = p->x;
    p_y = p->y;

    float spacing = params->grid_spacing;

    // Calculate coordinates within bucket grid
    grid_x = floor(p_x/spacing);
    grid_y = floor(p_y/spacing);

    // Go through neighboring buckets
    for(int dy=-1; dy<=1; dy++) {
        for(int dx=-1; dx<=1; dx++) {

            // If the neighbor bucket is outside of the grid we don't process it
            if ( grid_y+dy < 0 || grid_x+dx < 0 || (grid_x+dx) >= params->grid_size_x || (grid_y+dy) >= params->grid_size_y)
                continue;

             // Linear hash index for bucket
             bucket_index = (grid_y+dy) *params->grid_size_x + grid_x+dx;

             // Start index for hash value of current neighbor bucket
             start_index = start_indexes[bucket_index];

             // If neighbor bucket is not empty
             if (start_index != 0xffffffff)
             {
                end_index = end_indexes[bucket_index];

                for(int j=start_index; j<end_index; j++)
                {
                    q = fluid_particle_pointers[particle_ids[j]];

                    QmP_x = (q->x-p_x);
                    QmP_y = (q->y-p_y);
                    r = sqrt(QmP_x*QmP_x + QmP_y*QmP_y);

                    ratio = r*h_recip;

                    OmR2 = (1.0f-ratio)*(1.0f-ratio); // (one - r)^2

                    if(ratio < 1.0f) {
                        p->density += OmR2;
                        p->density_near += OmR2*(1.0f-ratio);
                    }
                }
            }
        }
    }
}

__device__ void boundaryConditions(fluid_particle *p, AABB_t *boundary, param *params)
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
        }

    }

    // Make sure object is not outside boundary
    // The particle must not be equal to boundary max or hash potentially won't pick it up
    // as the particle will in the 'next' after last bin
    if(p->x < boundary->min_x) {
        p->x = boundary->min_x;
    }
    else if(p->x > boundary->max_x){
        p->x = boundary->max_x-0.001f;
    }
    if(p->y <  boundary->min_y) {
        p->y = boundary->min_y;
    }
    else if(p->y > boundary->max_y){
        p->y = boundary->max_y-0.001f;
    }
}

__device__ void checkVelocity(float *v_x, float *v_y)
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

__device__ uint hash_val(float x, float y, param *params)
{
    float spacing = params->grid_spacing;
    float size_x  = params->grid_size_x;

    // Calculate grid coordinates
    uint grid_x,grid_y;
    grid_x = floor(x/spacing);
    grid_y = floor(y/spacing);

    uint grid_position = (grid_y * size_x + grid_x);

    return grid_position;
}

// The following kernel is modified from CUDA SDK particles example
// rearrange particle data into sorted order, and find the start of each cell
// in the sorted hash array
__global__ void find_cell_start(uint   *start_indexes,        // output: cell start index
                                  uint   *end_indexes,          // output: cell end index
                                  uint   *hash_values, // input: sorted grid hashes
                                  uint   *particle_ids,// input: sorted particle indices
                                  uint    numParticles)
{
    extern __shared__ uint sharedHash[];    // blockSize + 1 elements
    uint index = __umul24(blockIdx.x,blockDim.x) + threadIdx.x;

    uint hash;

    // handle case when no. of particles not multiple of block size
    if (index < numParticles)
    {
        hash = hash_values[index];

        // Load hash data into shared memory so that we can look
        // at neighboring particle's hash value without loading
        // two hash values per thread
        sharedHash[threadIdx.x+1] = hash;

        if (index > 0 && threadIdx.x == 0)
        {
            // first thread in block must load neighbor particle hash
            sharedHash[0] = hash_values[index-1];
        }
    }

    __syncthreads();

    if (index < numParticles)
    {
        // If this particle has a different cell index to the previous
        // particle then it must be the first particle in the cell,
        // so store the index of this particle in the cell.
        // As it isn't the first particle, it must also be the cell end of
        // the previous particle's cell

        if (index == 0 || hash != sharedHash[threadIdx.x])
        {
            start_indexes[hash] = index;

            if (index > 0)
                end_indexes[sharedHash[threadIdx.x]] = index;
        }

        if (index == numParticles - 1)
        {
            end_indexes[hash] = index + 1;
        }
    }

    // Potentially could allocate new_particles_array and reorder as in example
}

__global__ void calculate_hash(fluid_particle **fluid_particle_pointers, uint *hash_values, uint *particle_ids, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    fluid_particle *p;
    int num_fluid = params->number_fluid_particles_local + params->number_halo_particles;

    if(i >= num_fluid)
        return;

    p = fluid_particle_pointers[i];
    hash_values[i] =  hash_val(p->x, p->y, params);
    particle_ids[i] = i;
}

__global__ void apply_gravity(fluid_particle **fluid_particle_pointers, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;

    fluid_particle *p;
    float dt = params->tunable_params.time_step;
    float g = -params->tunable_params.g;

    int num_fluid = params->number_fluid_particles_local + params->number_halo_particles;

    if(i >= num_fluid)
        return;

    p = fluid_particle_pointers[i];
    p->v_y += g*dt;

    // Zero out density as well
    p->density = 0.0f;
    p->density_near = 0.0f;
}

__global__ void viscosity_impluses(fluid_particle **fluid_particle_pointers, uint *particle_ids, uint *start_indexes, uint *end_indexes, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;

    int num_fluid, grid_x, grid_y, bucket_index;
    uint start_index, end_index;
    fluid_particle *p, *q;
    float r, r_recip, ratio, u, imp, imp_x, imp_y;
    float p_x, p_y;
    float QmP_x, QmP_y;
    float h_recip, sigma, beta, dt,spacing;

    num_fluid = params->number_fluid_particles_local;
    h_recip = 1.0f/params->tunable_params.smoothing_radius;
    sigma = params->tunable_params.sigma;
    beta = params->tunable_params.beta;
    dt = params->tunable_params.time_step;
    spacing = params->grid_spacing;

    if(i >= num_fluid)
        return;

    p = fluid_particle_pointers[i];
    p_x = p->x;
    p_y = p->y;

    // Calculate coordinates within bucket grid
    grid_x = floor(p_x/spacing);
    grid_y = floor(p_y/spacing);

    // Go through neighboring buckets
    for(int dy=-1; dy<=1; dy++) {
        for(int dx=-1; dx<=1; dx++) {

            // If the neighbor bucket is outside of the grid we don't process it
            if ( grid_y+dy < 0 || grid_x+dx < 0 || (grid_x+dx) >= params->grid_size_x || (grid_y+dy) >= params->grid_size_y)
                continue;

             // Linear hash index for bucket
             bucket_index = (grid_y+dy) *params->grid_size_x + grid_x+dx;

             // Start index for hash value of current neighbor bucket
             start_index = start_indexes[bucket_index];

             // If neighbor bucket is not empty
             if (start_index != 0xffffffff)
             {
                end_index = end_indexes[bucket_index];

                for(int j=start_index; j<end_index; j++)
                {
                    q = fluid_particle_pointers[particle_ids[j]];

                    // Continue if same particle
                    if (p==q)
                        continue;

                    QmP_x = (q->x-p_x);
                    QmP_y = (q->y-p_y);
                    r = sqrt(QmP_x*QmP_x + QmP_y*QmP_y);

                    r_recip = 1.0f/r;
                    ratio = r*h_recip;

                    //Inward radial velocity
                    u = ((p->v_x-q->v_x)*QmP_x + (p->v_y-q->v_y)*QmP_y)*r_recip;

                    if(u>0.0f && u<=1.0f)
                    {
                        imp = dt * (1-ratio)*(sigma * u + beta * u*u);
                        imp_x = imp*QmP_x*r_recip;
                        imp_y = imp*QmP_y*r_recip;

                        // Not correct to use velocity check but will stop velocity from
                        // blowing up
                        checkVelocity(&imp_x, &imp_y);

//                        p->v_x -= imp_x*0.5f;
//                        p->v_y -= imp_y*0.5f;
                    }
                 } // End neighbor bucket particle loop  

             } // bucket not empty

        } // end x
    }  // end y

}

__global__ void predict_positions(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int num_fluid = params->number_fluid_particles_local;
    fluid_particle *p;
    float dt = params->tunable_params.time_step;

    if(i >= num_fluid)
        return;
    p = fluid_particle_pointers[i];
    p->x_prev = p->x;
    p->y_prev = p->y;
    p->x += (p->v_x * dt);
    p->y += (p->v_y * dt);

    // Enforce boundary conditions
    boundaryConditions(p, boundary_global, params);
}

__device__ void updateVelocity(fluid_particle *p, param *params)
{
    float dt = params->tunable_params.time_step;
    float v_x, v_y;

    v_x = (p->x-p->x_prev)/dt;
    v_y = (p->y-p->y_prev)/dt;

    checkVelocity(&v_x, &v_y);

    p->v_x = v_x;
    p->v_y = v_y;
}

__global__ void updateVelocities(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int num_fluid = params->number_fluid_particles_local;

    if(i >= num_fluid)
        return;

    fluid_particle *p;
    p = fluid_particle_pointers[i];
    boundaryConditions(p, boundary_global, params);
    updateVelocity(p, params);
}

__global__ void calculate_pressure(fluid_particle **fluid_particle_pointers, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int num_fluid = params->number_fluid_particles_local + params->number_halo_particles;

    float k, k_near, rest_density;

    k = params->tunable_params.k;
    k_near = params->tunable_params.k_near;
    rest_density = params->tunable_params.rest_density;

    fluid_particle *p;

    if(i >= num_fluid)
        return;

    p = fluid_particle_pointers[i];
    // Compute pressure and near pressure
    p->pressure = k * (p->density - rest_density);
    p->pressure_near = k_near * p->density_near;
}

__global__ void double_density_relaxation(fluid_particle **fluid_particle_pointers, uint *particle_ids, uint *start_indexes, uint *end_indexes, param *params)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int num_fluid = params->number_fluid_particles_local;

    if(i >= num_fluid)
        return;

    int bucket_index, start_index, end_index;
    fluid_particle *p, *q;
    float r,ratio,dt,h,h_recip,r_recip,D,D_x,D_y;
    float k_spring, p_pressure, p_pressure_near;
    float OmR, grid_x, grid_y;

    num_fluid = params->number_fluid_particles_local;
    k_spring = params->tunable_params.k_spring;
    h = params->tunable_params.smoothing_radius;
    h_recip = 1.0f/h;
    dt = params->tunable_params.time_step;

    // Iterating through the array in reverse reduces biased particle movement
    p = fluid_particle_pointers[i];
    p_pressure = p->pressure;
    p_pressure_near = p->pressure_near;

    float spacing = params->grid_spacing;

    // Calculate coordinates within bucket grid
    grid_x = floor(p->x/spacing);
    grid_y = floor(p->y/spacing);

    // Go through neighboring buckets
    for(int dy=-1; dy<=1; dy++) {
        for(int dx=-1; dx<=1; dx++) {

            // If the neighbor bucket is outside of the grid we don't process it
            if ( grid_y+dy < 0 || grid_x+dx < 0 || (grid_x+dx) >= params->grid_size_x || (grid_y+dy) >= params->grid_size_y)
                continue;

             // Linear hash index for bucket
             bucket_index = (grid_y+dy) *params->grid_size_x + grid_x+dx;

             // Start index for hash value of current neighbor bucket
             start_index = start_indexes[bucket_index];

             // If neighbor bucket is not empty
             if (start_index != 0xffffffff)
             {
                end_index = end_indexes[bucket_index];

                for(int j=start_index; j<end_index; j++)
                {
                    q = fluid_particle_pointers[particle_ids[j]];

                    // Continue if same particle
                    if (p==q)
                        continue;

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

                        p->x -= D_x;
                        p->y -= D_y;
                  } // If in ratio
              }
            }
           }
      } 
}

extern "C" void double_density_relaxation_gpu(fluid_particle **fluid_particle_pointers, uint *particle_ids, uint *start_indexes, uint *end_indexes, param *params)
{
    int total_particles = params->number_fluid_particles_local;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    double_density_relaxation<<<num_blocks, block_size>>>(fluid_particle_pointers, particle_ids, start_indexes, end_indexes, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}

extern "C" void calculate_pressures_gpu(fluid_particle **fluid_particle_pointers, param *params)
{
    int total_particles = params->number_fluid_particles_local + params->number_halo_particles;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    calculate_pressure<<<num_blocks, block_size>>>(fluid_particle_pointers, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );


}

extern "C" void updateVelocities_gpu(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params)
{
    int total_particles = params->number_fluid_particles_local;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    updateVelocities<<<num_blocks, block_size>>>(fluid_particle_pointers, boundary_global, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}

extern "C" void predict_positions_gpu(fluid_particle **fluid_particle_pointers, AABB_t *boundary_global, param *params)
{
    int total_particles = params->number_fluid_particles_local;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    predict_positions<<<num_blocks, block_size>>>(fluid_particle_pointers, boundary_global, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}

// Use thrust radix sort to sort ()
// Could also use uint2 and leap iterator...
extern "C" void sort_hash_gpu(uint *d_particle_ids, uint *d_hash_values, param *params)
{
    int total_particles = params->number_fluid_particles_local + params->number_halo_particles;

    thrust::device_ptr<uint> keys(d_hash_values);
    thrust::device_ptr<uint> values(d_particle_ids);

    thrust::sort_by_key(keys, keys+total_particles, values);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

extern "C" void hash_particles_gpu(fluid_particle **fluid_particle_pointers, uint *hash_values, uint *particle_ids, uint *starts, uint *ends, param *params)
{
    int total_particles = params->number_fluid_particles_local + params->number_halo_particles;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    // Reset start indexes
    unsigned int length_hash = params->grid_size_x * params->grid_size_y;
    cudaMemset(starts, 0xffffffff, length_hash*sizeof(uint));

    // Hash particles
    calculate_hash<<<num_blocks, block_size>>>(fluid_particle_pointers, hash_values, particle_ids, params);  

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

    // Sort hashed values
    sort_hash_gpu(particle_ids, hash_values, params);

    // Find start/end indexes for sorted values
    uint smem_size = sizeof(uint)*(block_size+1);
    num_blocks = ceil(length_hash/(float)block_size);
    find_cell_start<<<num_blocks, block_size, smem_size>>>(starts, ends, hash_values, particle_ids, total_particles);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}

extern "C" void apply_gravity_gpu(fluid_particle **fluid_particle_pointers, param *params)
{
    int total_particles = params->number_fluid_particles_local + params->number_halo_particles;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    apply_gravity<<< num_blocks, block_size >>>(fluid_particle_pointers, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );
}

extern "C" void calculate_density_gpu(fluid_particle **fluid_particle_pointers, uint *start_indexes, uint *end_indexes, uint *particle_ids, param *params)
{
    int total_particles = params->number_fluid_particles_local + params->number_halo_particles;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    calculate_density<<<num_blocks, block_size>>>(fluid_particle_pointers, start_indexes, end_indexes, particle_ids, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}

extern "C" void viscosity_impluses_gpu(fluid_particle **fluid_particle_pointers, uint *particle_ids, uint *start_indexes, uint *end_indexes, param *params)
{
    int total_particles = params->number_fluid_particles_local;
    int block_size = 256;
    int num_blocks = ceil(total_particles/(float)block_size);

    viscosity_impluses<<<num_blocks, block_size>>>(fluid_particle_pointers, particle_ids, start_indexes, end_indexes, params);

    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaDeviceSynchronize() );

}
