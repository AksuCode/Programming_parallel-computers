#include <math.h>
#include <vector>

#include <iostream>

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/
typedef double double4_t
__attribute__ ((vector_size (4 * sizeof(double))));

void correlate(int ny, int nx, const float *data, float *result) {

        int partition_count = 4;
        int partition_step = ceil(double(ny)/double(partition_count));

        for (int n = 0; n < ny && n < partition_count; n++) {
            int partition_limit = (n + 1) * partition_step;
            if (partition_limit >= ny) partition_limit = ny;
            for (int j = n * partition_step; j < partition_limit; j++) {

                int l = 0;
                for (int i = 0; i < nx; i++) {

                    l++;
                    if (l >= 4) l = 0;
                }


            }
        };


}

/*
    int core_count = 4;
    int thread_partition = floor(double(ny)/double(core_count));
    if (core_count > ny) thread_partition=1;

    for (int n = 0; n < ny && n < core_count; n++) {
        for (int i = n * thread_partition; i < ny && i < (n + 1) * thread_partition; i++) {
            for (int j = (n + 1) * thread_partition; j >= i; j--) {
                double4_t sum = {0, 0, 0, 0};
                for (int k = 0; k < vector_c_for_row; k++) {
                    double4_t prod = vectorized_m[k + j * vector_c_for_row] * vectorized_m[k + i * vector_c_for_row];
                    sum = sum + prod;
                }
                double final_sum = 0;
                for (int n = 0; n < 4; n++) {
                    final_sum = final_sum + sum[n];
                }
                result[i + ny * j] = final_sum;
            }
        }

    }
*/

/*
void correlate(int ny, int nx, const float *data, float *result) {

    int vector_c_for_row = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m = std::vector<double4_t> (ny * vector_c_for_row, double4_t {0, 0, 0, 0});
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        int l = 0;
        int n = 0;
        double row_mean = 0;
        for (int i = 0; i < nx; i++) {
            double new_dat_p = data[i + nx * j];
            vectorized_m[n + j * vector_c_for_row][l] = new_dat_p;
            row_mean = row_mean + new_dat_p;
            l++;
            if (l>=4) {
                l = 0;
                n++;
            }
        }

        row_mean = row_mean/nx;
        double4_t vectorized_row_mean = {row_mean, row_mean, row_mean, row_mean};
        double4_t padded_vectorized_row_mean = {row_mean, row_mean, row_mean, row_mean};
        if(l == 0) l = 4;
        for (; l < 4; l++) {
            padded_vectorized_row_mean[l] = 0;
        }
        int k = 0;
        for (; k < vector_c_for_row - 1; k++) {
            vectorized_m[k + j * vector_c_for_row] = vectorized_m[k + j * vector_c_for_row] - vectorized_row_mean;
        }
        vectorized_m[k + j * vector_c_for_row] = vectorized_m[k + j * vector_c_for_row] - padded_vectorized_row_mean;
    }

    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        double4_t vectorized_sqr_sum = {0,0,0,0};
        for (int k = 0; k < vector_c_for_row; k++) {
            double4_t sqr = vectorized_m[k + j * vector_c_for_row] * vectorized_m[k + j * vector_c_for_row];
            vectorized_sqr_sum = vectorized_sqr_sum + sqr;
        }

        double sqr_sum = 0;
        for (int t = 0; t < 4; t++) {
            sqr_sum = sqr_sum + vectorized_sqr_sum[t];
        }

        double sqrt_sum = sqrt(sqr_sum);
        for (int k = 0; k < vector_c_for_row; k++) {
            vectorized_m[k + j * vector_c_for_row] = vectorized_m[k + j * vector_c_for_row]/sqrt_sum;
        }
    }

    int core_count = 4;
    int thread_partition = ceil(double(ny)/double(core_count));
    if (core_count > ny) thread_partition=1;

    for (int n = 0; n < ny && n < core_count; n++) {
        for (int i = n * thread_partition; i < ny && i < (n + 1) * thread_partition; i++) {


            int j = (n + 1) * thread_partition;
            if (j >= ny) j = ny - 1;
            for (; j >= i; j--) {
                double4_t sum = {0, 0, 0, 0};
                for (int k = 0; k < vector_c_for_row; k++) {
                    double4_t prod = vectorized_m[k + j * vector_c_for_row] * vectorized_m[k + i * vector_c_for_row];
                    sum = sum + prod;
                }
                double final_sum = 0;
                for (int n = 0; n < 4; n++) {
                    final_sum = final_sum + sum[n];
                }
                result[i + ny * j] = final_sum;
            }
        }

    }

}
*/