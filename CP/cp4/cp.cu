#include <vector>
#include <math.h>

#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>

static inline void check(cudaError_t err, const char* context) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA error: " << context << ": "
            << cudaGetErrorString(err) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

#define CHECK(x) check(x, #x)

static inline int divup(int a, int b) {
    return (a + b - 1)/b;
}

/*
static inline int roundup(int a, int b) {
    return divup(a, b) * b;
}
*/

__global__ void mykernel(float* result, float * normalized_data, int nx, int ny) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if (i >= ny || j >= ny)
        return;
    if (j > i) {
        result[i + ny * j] = 0;
        return;
    }
    float sum = 0;
    for (int k = 0; k < nx; k++) {
        sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
    }
    result[i + ny * j] = sum;
}

void correlate(int ny, int nx, const float *data, float *result) {

    std::vector<float> normalized_data = std::vector<float> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean/nx;
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = data[i + row_indx] - mean;
        }
        float sqr_sum = 0;
        for (int i = 0; i < nx; i++) {
            sqr_sum = sqr_sum + normalized_data[i + row_indx] * normalized_data[i + row_indx];
        }
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = normalized_data[i + row_indx]/sqrt(sqr_sum);
        }
    }

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, nx * ny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));
    CHECK(cudaMemcpy(dGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    dim3 dimBlock(16, 16);
    dim3 dimGrid(divup(ny, dimBlock.x), divup(ny, dimBlock.y));
    mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, ny);
    CHECK(cudaGetLastError());

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}