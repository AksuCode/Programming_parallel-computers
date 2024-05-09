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

    int vector_c_for_row = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m = std::vector<double4_t> (ny * vector_c_for_row, double4_t {0, 0, 0, 0});
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

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
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

/*
void correlate(int ny, int nx, const float *data, float *result) {

    int vector_c_for_row = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m = std::vector<double4_t> (ny * vector_c_for_row, double4_t {0, 0, 0, 0});
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

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
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

    std::vector<double> row_means = std::vector<double> (ny, 0);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            row_means[j] = row_means[j] + data[i + nx * j];
        }
        row_means[j] = row_means[j]/nx;
    }

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            normalized_data[i + nx * j] = data[i + nx * j] - row_means[j];
        }
    }

    std::vector<double> sqrd_sums = std::vector<double> (ny, 0);
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            sqrd_sums[j] = sqrd_sums[j] + normalized_data[i + nx * j] * normalized_data[i + nx * j];
        }
    }

    int vector_c_for_row = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_normal_m = std::vector<double4_t> (ny * vector_c_for_row, double4_t {0, 0, 0, 0});
    for (int j = 0; j < ny; j++) {
        int l = 0;
        int k = 0;
        for (int i = 0; i < nx; i++) {
            vectorized_normal_m[k + j * vector_c_for_row][l] = normalized_data[i + nx * j]/sqrt(sqrd_sums[j]);
            l++;
            if (l>=4) {
                l = 0;
                k++;
            }
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double4_t sum = {0, 0, 0, 0};
            for (int k = 0; k < vector_c_for_row; k++) {
                double4_t prod = vectorized_normal_m[k + j * vector_c_for_row] * vectorized_normal_m[k + i * vector_c_for_row];
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