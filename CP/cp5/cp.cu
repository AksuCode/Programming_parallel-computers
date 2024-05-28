#include <vector>
#include <math.h>

#include <cstdlib>
#include <iostream>
#include <cuda_runtime.h>

#include <iomanip>

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

static inline int roundup(int a, int b) {
    return divup(a, b) * b;
}

__global__ void mykernel(float* r, const float * d, int nx, int ny) {

    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x * 512;
    int jc = blockIdx.y * 512;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }

    for (int k = 0; k < nx; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic + ib * 8 + ia;
            if (i < ny) x[ib] = d[k + nx * i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc + jb * 8 + ja;
            if (j < ny) y[jb] = d[k + nx * j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic + ib * 8 + ia;
            int j = jc + jb * 8 + ja;
            if (i < ny && j < ny) {
                r[i + ny * j] = v[ib][jb];
            }
        }
    }

}

__global__ void myppkernel(const float* r, float* d, int nx, int ny) {

    int j = threadIdx.x + blockIdx.x * blockDim.x;

    double divisor = 1/double(nx);

    int row_indx = nx * j;
    double mean = 0;
    if (j < ny) {
        for (int i = 0; i < nx; i++) {
            mean = mean + double(r[i + row_indx]);
        }
        mean = mean * divisor;
        for (int i = 0; i < nx; i++) {
            d[i + row_indx] = float(double(r[i + row_indx]) - mean);
        }
        double sqr_sum = 0;
        for (int i = 0; i < nx; i++) {
            sqr_sum = sqr_sum + double(d[i + row_indx]) * double(d[i + row_indx]);
        }
        for (int i = 0; i < nx; i++) {
            d[i + row_indx] = double(d[i + row_indx])/double(sqrt(sqr_sum));
        }
    }

}

void correlate(int ny, int nx, const float *data, float *result) {

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, nx * ny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, nx * ny * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, data, nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64);
        dim3 dimGrid(roundup(ny, dimBlock.x));
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, ny);
        CHECK(cudaGetLastError());
    }

    CHECK(cudaFree(rGPU));
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(roundup(ny, 512), roundup(ny, 512));
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, ny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}







/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }

    for (int k = 0; k < nx; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic * 64 + ib * 8 + ia;
            x[ib] = d[k + nnx * i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc * 64 + jb * 8 + ja;
            y[jb] = d[k + nnx * j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            if (i < ny && j < ny) {
                r[nx*i + j] = v[ib][jb];
            }
        }
    }
}

__global__ void myppkernel(const float* r, float* d, int nx, int nnx, int ny) {
    int ja = threadIdx.x;
    int i = blockIdx.y;

    for (int jb = 0; jb < nnx; jb += 64) {
        int j = jb + ja;
        float v = (i < nx && j < ny) ? r[nx*i + j] : 0;
        d[nnx*i + j] = v;
    }
}

void correlate(int ny, int nx, const float *data, float *result) {

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    std::vector<float> normalized_data(nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int k = 0; k < nx; k++) {
            mean = mean + data[k + row_indx];
        }
        mean = mean/nx;
        for (int k = 0; k < nx; k++) {
            normalized_data[k + row_indx] = data[k + row_indx] - mean;
        }
        float sqr_sum = 0;
        for (int k = 0; k < nx; k++) {
            sqr_sum = sqr_sum + normalized_data[k + row_indx] * normalized_data[k + row_indx];
        }
        for (int k = 0; k < nx; k++) {
            normalized_data[k + row_indx] = normalized_data[k + row_indx]/sqrt(sqr_sum);
        }
    }

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, nx * ny * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nnx);
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny);
        CHECK(cudaGetLastError());
    }

    CHECK(cudaFree(rGPU));
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nny / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));
}
*/





/*
__global__ void mykernel(float* r, const float* d, int n, int nn) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }

    for (int k = 0; k < n; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic * 64 + ib * 8 + ia;
            x[ib] = d[k + nn * i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc * 64 + jb * 8 + ja;
            y[jb] = d[k + nn * j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            if (i < n && j < n) {
                r[n*i + j] = v[ib][jb];
            }
        }
    }
}

__global__ void myppkernel(const float* r, float* d, int n, int nn) {
    int ja = threadIdx.x;
    int i = blockIdx.y;

    for (int jb = 0; jb < nn; jb += 64) {
        int j = jb + ja;
        float v = (i < n && j < n) ? r[n*i + j] : 0;
        d[nn*i + j] = v;
    }
}

void correlate(int ny, int nx, const float *data, float *result) {

    int n = nx;
    int nn = roundup(n, 64);

    std::vector<float> normalized_data(nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = n * j;
        float mean = 0;
        for (int k = 0; k < n; k++) {
            mean = mean + data[k + row_indx];
        }
        mean = mean/n;
        for (int k = 0; k < n; k++) {
            normalized_data[k + row_indx] = data[k + row_indx] - mean;
        }
        float sqr_sum = 0;
        for (int k = 0; k < n; k++) {
            sqr_sum = sqr_sum + normalized_data[k + row_indx] * normalized_data[k + row_indx];
        }
        for (int k = 0; k < n; k++) {
            normalized_data[k + row_indx] = normalized_data[k + row_indx]/sqrt(sqr_sum);
        }
    }

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, nn * nn * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, n * n * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, &normalized_data[0], n * n * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nn);
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, n, nn);
        CHECK(cudaGetLastError());
    }

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nn / 64, nn / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, n, nn);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, n * n * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));
}
*/



/*
__global__ void mykernel(float* r, const float* d, int n, int nn) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    const float* t = d + nn * nn;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }

    for (int k = 0; k < n; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic * 64 + ib * 8 + ia;
            x[ib] = t[nn*k + i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc * 64 + jb * 8 + ja;
            y[jb] = d[nn*k + j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            if (i < n && j < n) {
                r[n*i + j] = v[ib][jb];
            }
        }
    }
}

__global__ void myppkernel(const float* r, float* d, int n, int nn) {
    int ja = threadIdx.x;
    int i = blockIdx.y;

    float* t = d + nn * nn;

    for (int jb = 0; jb < nn; jb += 64) {
        int j = jb + ja;
        float v = (i < n && j < n) ? r[n*i + j] : 0;
        d[nn*i + j] = v;
        t[nn*j + i] = v;
    }
}

void correlate(int ny, int nx, const float *data, float *result) {

    int n = nx;
    int nn = roundup(n, 64);

    std::vector<float> normalized_data(nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = n * j;
        float mean = 0;
        for (int k = 0; k < n; k++) {
            mean = mean + data[k + row_indx];
        }
        mean = mean/n;
        for (int k = 0; k < n; k++) {
            normalized_data[k + row_indx] = data[k + row_indx] - mean;
        }
        float sqr_sum = 0;
        for (int k = 0; k < n; k++) {
            sqr_sum = sqr_sum + normalized_data[k + row_indx] * normalized_data[k + row_indx];
        }
        for (int k = 0; k < n; k++) {
            normalized_data[k + row_indx] = normalized_data[k + row_indx]/sqrt(sqr_sum);
        }
    }

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, 2 * nn * nn * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, n * n * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, &normalized_data[0], n * n * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nn);
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, n, nn);
        CHECK(cudaGetLastError());
    }

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nn / 64, nn / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, n, nn);
        CHECK(cudaGetLastError());
    }

    float sus[2 * nn * nn];
    CHECK(cudaMemcpy(sus, dGPU, 2 * nn * nn * sizeof(float), cudaMemcpyDeviceToHost));

    std::cout << "Padded Ndata result: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << sus[i + nn * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    const float* t = sus + nn * nn;

    std::vector<float> sus_res(nx * ny, 0);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; ++k) {
                float x;
                float y;
                x = t[j + nn * k];
                y = sus[k + nn * i];
                sus_res[n*i + j] = sus_res[n*i + j] + x * y;
            }
        }
    }

    std::cout << "Result from kernel: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << sus_res[i + n * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Ndata result: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << normalized_data[i + nx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    float linear_res[nx * ny];
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < ny; j++) {
            double sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            linear_res[i + ny * j] = sum;
        }
    }
    std::cout << "Linear result: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << linear_res[i + nx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;


    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, n * n * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));
}
*/
















































/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    const float* t = d + nnx * nny;

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            float sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + d[k + nnx * i] * t[j + nnx * k];
            }
            if (i < ny && j < ny) {
                if (sum == 0) r[ny*i + j] = 69;
                else r[ny*i + j] = sum;
            }
        }
    }

}
*/

/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    const float* t = d + nnx * nny;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }

    for (int k = 0; k < nx; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic * 64 + ib * 8 + ia;
            x[ib] = t[nnx*k + i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc * 64 + jb * 8 + ja;
            y[jb] = d[nnx*k + j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            if (i < ny && j < ny) {
                r[ny*i + j] = v[ib][jb];
            }
        }
    }

}


__global__ void myppkernel(const float* r, float* d, int nx, int nnx, int ny, int nny) {
    int ja = threadIdx.x;
    int i = blockIdx.y;

    float* t = d + nnx * nny;

    for (int jb = 0; jb < nnx; jb += 64) {
        int j = jb + ja;
        float v = (i < nx && j < ny) ? r[nx*i + j] : 0;
        d[nnx*i + j] = v;
        t[nnx*j + i] = v;
    }
}

void correlate(int ny, int nx, const float *data, float *result) {

    std::vector<float> normalized_data = std::vector<float> (nx * ny);
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

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, 2 * nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, nx * ny * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nnx);   //?
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    CHECK(cudaFree(rGPU));
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nnx / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));
}
*/



/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    const float* t = d + nnx * nny;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }

    for (int k = 0; k < nx; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic * 64 + ib * 8 + ia;
            x[ib] = t[nnx*k + i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc * 64 + jb * 8 + ja;
            y[jb] = d[nnx*k + j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            if (i < ny && j < ny) {
                r[ny*i + j] = v[ib][jb];
            }
        }
    }
}

__global__ void myppkernel(const float* r, float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int j = blockIdx.y;

    float* t = d + nnx * nny;

    for (int ib = 0; ib < nnx; ib += 64) {
        int i = ib + ia;
        float v = (i < nx && j < ny) ? r[i + nx * j] : 0;
        d[i + nnx * j] = v;
        t[j + nnx * i] = v;
    }

    __syncthreads();

    float mean = 0;
    for (int i = 0; i < nx; i++) {
        mean = mean + d[i + nnx * j];
    }
    mean = mean/nx;

    for (int i = 0; i < nx; i++) {
        d[i + nnx * j] = d[i + nnx * j] - mean;
        t[j + nnx * i] = t[j + nnx * i] - mean;
    }

    float sqr_sum = 0;
    for (int i = 0; i < nx; i++) {
        sqr_sum = sqr_sum + d[i + nnx * j] * d[i + nnx * j];
    }

    for (int i = 0; i < nx; i++) {
        d[i + nnx * j] = d[i + nnx * j]/sqrt(sqr_sum);
        d[j + nnx * i] = d[j + nnx * i]/sqrt(sqr_sum);
    }

}


void correlate(int ny, int nx, const float *data, float *result) {

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, 2 * nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, nx * ny * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, data, nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nny);
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }


    float sus [2 * nnx * nny];
    float* t = sus + nnx * nny;
    CHECK(cudaMemcpy(sus, dGPU, 2 * nnx * nny * sizeof(float), cudaMemcpyDeviceToHost));
    std::cout << "GPU ndata: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << sus[i + nnx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "GPU T ndata: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << t[i + nnx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;


    float divisor = 1/float(nx);

    std::vector<float> normalized_data = std::vector<float> (nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean * divisor;
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


    std::cout << "Linear ndata: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << normalized_data[i + nx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;


    CHECK(cudaFree(rGPU));
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nnx / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}
*/





/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            float sum = 0;
            for (int k = 0; k < nnx; k++) {
                sum = sum + d[k + nnx * i] * d[k + nnx * j];
            }
            if (i < ny && j < ny) {
                r[ny*i + j] = sum;
            }
        }
    }
}


__global__ void myppkernel(float* d, const float* data, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int j = blockIdx.y;

    for (int ib = 0; ib < nnx; ib += 64) {
        int i = ia + ib;
        float v = (i < nx && j < ny) ? data[i + nx * j] : 0;
        d[i + nnx * j] = v;
    }
}


void correlate(int ny, int nx, const float *data, float *result) {

    float divisor = 1/float(nx);

    std::vector<float> normalized_data = std::vector<float> (nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean * divisor;
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

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    // Allocate memory & copy data to GPU
    float* dataGPU = NULL;
    CHECK(cudaMalloc((void**)&dataGPU, nx * ny * sizeof(float)));
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));
    CHECK(cudaMemcpy(dataGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nnx);
        myppkernel<<<dimGrid, dimBlock>>>(dGPU, dataGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    CHECK(cudaFree(dataGPU));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nnx / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}
*/


/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    const float* t = d + nnx * nny;

    float v[8][8];
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            v[ib][jb] = 0;
        }
    }
    for (int k = 0; k < nx; ++k) {
        float x[8];
        float y[8];
        for (int ib = 0; ib < 8; ++ib) {
            int i = ic * 64 + ib * 8 + ia;
            x[ib] = t[nnx*k + i];
        }
        for (int jb = 0; jb < 8; ++jb) {
            int j = jc * 64 + jb * 8 + ja;
            y[jb] = d[nnx*k + j];
        }
        for (int ib = 0; ib < 8; ++ib) {
            for (int jb = 0; jb < 8; ++jb) {
                v[ib][jb] = v[ib][jb] + x[ib] * y[jb];
            }
        }
    }
    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            if (i < nx && j < ny) {
                r[ny*i + j] = v[ib][jb];
            }
        }
    }
}


__global__ void myppkernel(float* d, const float* r, int nx, int nnx, int ny, int nny) {
    int ja = threadIdx.x;
    int i = blockIdx.y;

    float* t = d + nnx * nny;

    for (int jb = 0; jb < nnx; jb += 64) {
        int j = jb + ja;
        float v = (i < nx && j < ny) ? r[nx*i + j] : 0;
        d[nnx*i + j] = v;
        t[nnx*j + i] = v;
    }
}


void correlate(int ny, int nx, const float *data, float *result) {

    float divisor = 1/float(nx);

    std::vector<float> normalized_data = std::vector<float> (nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean * divisor;
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

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    // Allocate memory & copy data to GPU
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, 2 * nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, nx * ny * sizeof(float)));
    CHECK(cudaMemcpy(rGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nnx);
        myppkernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    CHECK(cudaFree(rGPU));
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nnx / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}
*/


/*
__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            float sum = 0;
            for (int k = 0; k < nnx; k++) {
                sum = sum + d[k + nnx * i] * d[k + nnx * j];
            }
            if (i < nx && j < ny) {
                r[ny*i + j] = sum;
            }
        }
    }
}


__global__ void myppkernel(float* d, const float* data, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int j = blockIdx.y;

    for (int ib = 0; ib < nnx; ib += 64) {
        int i = ia + ib;
        float v = (i < nx && j < ny) ? data[i + nx * j] : 0;
        d[i + nnx * j] = v;
    }
}


void correlate(int ny, int nx, const float *data, float *result) {

    float divisor = 1/float(nx);

    std::vector<float> normalized_data = std::vector<float> (nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean * divisor;
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

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    // Allocate memory & copy data to GPU
    float* dataGPU = NULL;
    CHECK(cudaMalloc((void**)&dataGPU, nx * ny * sizeof(float)));
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));
    CHECK(cudaMemcpy(dataGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nnx);
        myppkernel<<<dimGrid, dimBlock>>>(dGPU, dataGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    float sus [nnx * ny];
    CHECK(cudaMemcpy(sus, dGPU, nnx * nny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dataGPU));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nnx / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}
*/


/*
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

static inline int roundup(int a, int b) {
    return divup(a, b) * b;
}

__global__ void mykernel(float* r, const float* d, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int ja = threadIdx.y;
    int ic = blockIdx.x;
    int jc = blockIdx.y;

    for (int ib = 0; ib < 8; ++ib) {
        for (int jb = 0; jb < 8; ++jb) {
            int i = ic * 64 + ib * 8 + ia;
            int j = jc * 64 + jb * 8 + ja;
            float sum = 0;
            for (int k = 0; k < nnx; k++) {
                sum = sum + d[k + nnx * i] * d[k + nnx * j];
            }
            if (i < nx && j < ny) {
                r[ny*i + j] = sum;
            }
        }
    }
}


__global__ void myppkernel(float* d, const float* data, int nx, int nnx, int ny, int nny) {
    int ia = threadIdx.x;
    int j = blockIdx.y;

    for (int ib = 0; ib < nnx; ib += 64) {
        int i = ia + ib;
        float v = (i < nx && j < ny) ? data[i + nx * j] : 0;
        d[i + nnx * j] = v;
    }
}


void correlate(int ny, int nx, const float *data, float *result) {

    float divisor = 1/float(nx);

    std::vector<float> normalized_data = std::vector<float> (nx * ny);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        float mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean * divisor;
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

    int nnx = roundup(nx, 64);
    int nny = roundup(ny, 64);

    // Allocate memory & copy data to GPU
    float* dataGPU = NULL;
    CHECK(cudaMalloc((void**)&dataGPU, nx * ny * sizeof(float)));
    float* dGPU = NULL;
    CHECK(cudaMalloc((void**)&dGPU, 2 * nnx * nny * sizeof(float)));
    float* rGPU = NULL;
    CHECK(cudaMalloc((void**)&rGPU, ny * ny * sizeof(float)));
    CHECK(cudaMemcpy(dataGPU, &normalized_data[0], nx * ny * sizeof(float), cudaMemcpyHostToDevice));

    // Run kernel
    {
        dim3 dimBlock(64, 1);
        dim3 dimGrid(1, nnx);
        myppkernel<<<dimGrid, dimBlock>>>(dGPU, dataGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }

    float sus [nnx * ny];
    CHECK(cudaMemcpy(sus, dGPU, nnx * ny * sizeof(float), cudaMemcpyDeviceToHost));
    std::cout << "GPU ndata: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << sus[i + nnx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Linear ndata: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << normalized_data[i + nx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    CHECK(cudaFree(dataGPU));

    // Run kernel
    {
        dim3 dimBlock(8, 8);
        dim3 dimGrid(nnx / 64, nny / 64);
        mykernel<<<dimGrid, dimBlock>>>(rGPU, dGPU, nx, nnx, ny, nny);
        CHECK(cudaGetLastError());
    }



    float linear_res[nx * ny];
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j < ny; j++) {
            double sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            linear_res[i + ny * j] = sum;
        }
    }
    std::cout << "Linear result: " << std::endl;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            std::cout.width(10);
            std::cout << std::fixed << std::setprecision(3) << linear_res[i + nx * j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Copy data back to CPU & release memory
    CHECK(cudaMemcpy(result, rGPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost));
    CHECK(cudaFree(dGPU));
    CHECK(cudaFree(rGPU));

}
*/