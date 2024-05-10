#include <math.h>
#include <vector>

#include <stdlib.h>
#include <chrono>
#include <algorithm>

#include <iostream>

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/

typedef double double8_t
__attribute__ ((vector_size (8 * sizeof(double))));

void correlate(int ny, int nx, const float *data, float *result) {

    // Ready constants:
    double dividend_by_size = 1/double(nx);
    int row_v_partition_count = ceil(double(nx)/double(8));

    std::vector<double8_t> vectorized_m = std::vector<double8_t> (ny * row_v_partition_count);
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {

        int n = 0;
        int l = 0;
        double8_t tmp_sum_v = double8_t {0,0,0,0,0,0,0,0};
        double8_t sum_v = double8_t {0,0,0,0,0,0,0,0};
        double8_t tmp_sqr_sum_v = double8_t {0,0,0,0,0,0,0,0};
        double8_t sqr_sum_v = double8_t {0,0,0,0,0,0,0,0};
        std::vector<double8_t> vectorized_row = std::vector<double8_t> (row_v_partition_count, double8_t {0,0,0,0,0,0,0,0});
        for (int i = 0; i < nx; i++) {
            double new_dat_p = data[i + nx * j];
            vectorized_row[n][l] = new_dat_p;
            tmp_sum_v[l] = new_dat_p;
            tmp_sqr_sum_v[l] = new_dat_p * new_dat_p;
            l++;
            if (l>=8) {
                sum_v = sum_v + tmp_sum_v;
                sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
                tmp_sum_v = double8_t {0,0,0,0,0,0,0,0};         // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                tmp_sqr_sum_v = double8_t {0,0,0,0,0,0,0,0};     // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                l = 0;
                n++;
            }
        }
        if(l != 0) {
            sum_v = sum_v + tmp_sum_v;
            sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
        }

        double sum_1 = sum_v[0] + sum_v[1];
        double sum_2 = sum_v[2] + sum_v[3];
        double sum_3 = sum_v[4] + sum_v[5];
        double sum_4 = sum_v[6] + sum_v[7];
        double sqr_sum_1 = sqr_sum_v[0] + sqr_sum_v[1];
        double sqr_sum_2 = sqr_sum_v[2] + sqr_sum_v[3];
        double sqr_sum_3 = sqr_sum_v[4] + sqr_sum_v[5];
        double sqr_sum_4 = sqr_sum_v[6] + sqr_sum_v[7];
        double sum_12 = sum_1 + sum_2;
        double sum_34 = sum_3 + sum_4;
        double sqr_sum_12 = sqr_sum_1 + sqr_sum_2;
        double sqr_sum_34 = sqr_sum_3 + sqr_sum_4;
        double sum = sum_12 + sum_34;
        double sqr_sum = sqr_sum_12 + sqr_sum_34;

        double mean = sum * dividend_by_size;

        double second_term = -2 * mean;
        double third_term = mean*mean;
        second_term = second_term * sum;
        third_term = third_term * nx;
        sqr_sum = sqr_sum + second_term;
        sqr_sum = sqr_sum + third_term;
        sqr_sum = sqrt(sqr_sum);
        sqr_sum  = 1/sqr_sum;

        double8_t norm_sqrt_sum_v = {sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum};
        double8_t mean_v = {mean, mean, mean, mean, mean, mean, mean, mean};
        double8_t padded_mean_v = {mean, mean, mean, mean, mean, mean, mean, mean};
        if(l == 0) l = 8;
        for (; l < 8; l++) {
            padded_mean_v[l] = 0;
        }
        int s = 0;
        for(; s < row_v_partition_count - 1; s++) {
            double8_t temp = vectorized_row[s]-mean_v;
            vectorized_row[s] = temp * norm_sqrt_sum_v;
        }
        double8_t temp = vectorized_row[s]-padded_mean_v;
        vectorized_row[s] = temp * norm_sqrt_sum_v;

        std::copy(vectorized_row.begin(), vectorized_row.end(), vectorized_m.begin() + (j * row_v_partition_count));
    }



    int block_size = 2;
    // Padding:
    int padding_size = ny % block_size;
    std::vector<double8_t> padding = std::vector<double8_t> (padding_size * row_v_partition_count, double8_t {0,0,0,0,0,0,0,0});
    vectorized_m.insert(vectorized_m.end(), padding.begin(), padding.end());
    int b_ny = ny + padding_size;
    int b_lim = b_ny/block_size;

    #pragma omp parallel for
    for (int ib = 0; ib < b_lim; ib++) {
        for (int jb = 0; jb < b_lim; jb++) {
            double8_t sum_v_00 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_10 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_01 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_11 = {0, 0, 0, 0, 0, 0, 0, 0};
            for (int k = 0; k < row_v_partition_count; k++) {
                double8_t prod_00 = vectorized_m[k + jb * row_v_partition_count] * vectorized_m[k + ib * row_v_partition_count];
                double8_t prod_10 = vectorized_m[k + jb * row_v_partition_count] * vectorized_m[k + (ib + 1) * row_v_partition_count];
                double8_t prod_01 = vectorized_m[k + (jb + 1) * row_v_partition_count] * vectorized_m[k + ib * row_v_partition_count];
                double8_t prod_11 = vectorized_m[k + (jb + 1) * row_v_partition_count] * vectorized_m[k + (ib + 1) * row_v_partition_count];
                sum_v_00 = sum_v_00 + prod_00;
                sum_v_10 = sum_v_10 + prod_10;
                sum_v_01 = sum_v_01 + prod_01;
                sum_v_11 = sum_v_11 + prod_11;
            }

            double sum_00 = 0;
            double sum_10 = 0;
            double sum_01 = 0;
            double sum_11 = 0;
            for (int k = 0; k < 8; k++) {
                sum_00 = sum_00 + sum_v_00[k];
                sum_10 = sum_10 + sum_v_10[k];
                sum_01 = sum_01 + sum_v_01[k];
                sum_11 = sum_11 + sum_v_11[k];
            }

            result[ib + ny * jb] = sum_00;
            result[(ib + 1) + ny * jb] = sum_10;
            result[ib + ny * (jb + 1)] = sum_01;
            result[(ib + 1) + ny * (jb + 1)] = sum_11;
        }

    }

}


/*
typedef double double4_t
__attribute__ ((vector_size (4 * sizeof(double))));

void correlate(int ny, int nx, const float *data, float *result) {

    auto t1 = std::chrono::high_resolution_clock::now();
    // Ready constants:
    double dividend_by_size = 1/double(nx);

    int row_v_partition_count = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m = std::vector<double4_t> (ny * row_v_partition_count);
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {

        int n = 0;
        int l = 0;
        double4_t tmp_sum_v = double4_t {0,0,0,0};
        double4_t sum_v = double4_t {0,0,0,0};
        double4_t tmp_sqr_sum_v = double4_t {0,0,0,0};
        double4_t sqr_sum_v = double4_t {0,0,0,0};
        std::vector<double4_t> vectorized_row = std::vector<double4_t> (row_v_partition_count, double4_t {0, 0, 0, 0});
        for (int i = 0; i < nx; i++) {
            double new_dat_p = data[i + nx * j];
            vectorized_row[n][l] = new_dat_p;
            tmp_sum_v[l] = new_dat_p;
            tmp_sqr_sum_v[l] = new_dat_p * new_dat_p;
            l++;
            if (l>=4) {
                sum_v = sum_v + tmp_sum_v;
                sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
                tmp_sum_v = double4_t{0,0,0,0};         // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                tmp_sqr_sum_v = double4_t{0,0,0,0};     // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                l = 0;
                n++;
            }
        }
        if(l != 0) {
            sum_v = sum_v + tmp_sum_v;
            sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
        }

        double sum = sum_v[0] + sum_v[1] + sum_v[2] + sum_v[3];
        double sqr_sum = sqr_sum_v[0] + sqr_sum_v[1] + sqr_sum_v[2] + sqr_sum_v[3];
        double mean = sum * dividend_by_size;

        double second_term = -2 * mean;
        double third_term = mean*mean;
        second_term = second_term * sum;
        third_term = third_term * nx;
        sqr_sum = sqr_sum + second_term;
        sqr_sum = sqr_sum + third_term;
        sqr_sum = sqrt(sqr_sum);
        sqr_sum  = 1/sqr_sum;

        double4_t norm_sqrt_sum_v = {sqr_sum, sqr_sum, sqr_sum, sqr_sum};
        double4_t mean_v = {mean, mean, mean, mean};
        double4_t padded_mean_v = {mean, mean, mean, mean};
        if(l == 0) l = 4;
        for (; l < 4; l++) {
            padded_mean_v[l] = 0;
        }
        int s = 0;
        for(; s < row_v_partition_count - 1; s++) {
            double4_t temp = vectorized_row[s]-mean_v;
            vectorized_row[s] = temp * norm_sqrt_sum_v;
        }
        double4_t temp = vectorized_row[s]-padded_mean_v;
        vectorized_row[s] = temp * norm_sqrt_sum_v;

        std::copy(vectorized_row.begin(), vectorized_row.end(), vectorized_m.begin() + (j * row_v_partition_count));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto ns_int_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    std::cout << ns_int_1.count() << std::endl;


    auto t3 = std::chrono::high_resolution_clock::now();
    int halved_row_v_partition_count = floor(double(row_v_partition_count)/double(2));
    #pragma omp parallel for schedule(static,1)
    for (int sus = 0; sus < ny; sus = sus + 24) {
        for (int i = 0; i < ny; i++) {
            for (int j = sus; j < sus + 24 && j <= i; j++) {
                double4_t sum = {0, 0, 0, 0};
                int k = 0;
                for (; k < halved_row_v_partition_count - 1; k++) {
                    double4_t prod_1 = vectorized_m[k*2 + j * row_v_partition_count] * vectorized_m[k*2 + i * row_v_partition_count];
                    double4_t prod_2 = vectorized_m[1 + k*2 + j * row_v_partition_count] * vectorized_m[1 + k*2 + i * row_v_partition_count];
                    double4_t tmp_sum = prod_1 + prod_2;
                    sum = sum + tmp_sum;
                }
                k = k * 2;
                for (; k < row_v_partition_count; k++) {
                    double4_t prod = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
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
    auto t4 = std::chrono::high_resolution_clock::now();
    auto ns_int_2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);
    std::cout << ns_int_2.count() << std::endl;


}
*/

/*
void correlate(int ny, int nx, const float *data, float *result) {

    int row_v_partition_count = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m = std::vector<double4_t> {};
    for (int j = 0; j < ny; j++) {

        int n = 0;
        int l = 0;
        double4_t tmp_sum_v = double4_t {0,0,0,0};
        double4_t sum_v = double4_t {0,0,0,0};
        double4_t tmp_sqr_sum_v = double4_t {0,0,0,0};
        double4_t sqr_sum_v = double4_t {0,0,0,0};
        std::vector<double4_t> vectorized_row = std::vector<double4_t> (row_v_partition_count, double4_t {0, 0, 0, 0});
        for (int i = 0; i < nx; i++) {
            double new_dat_p = data[i + nx * j];
            vectorized_row[n][l] = new_dat_p;
            tmp_sum_v[l] = new_dat_p;
            tmp_sqr_sum_v[l] = new_dat_p * new_dat_p;
            l++;
            if (l>=4) {
                sum_v = sum_v + tmp_sum_v;
                sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
                tmp_sum_v = double4_t{0,0,0,0};         // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                tmp_sqr_sum_v = double4_t{0,0,0,0};     // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                l = 0;
                n++;
            }
        }
        if(l != 0) {
            sum_v = sum_v + tmp_sum_v;
            sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
        }

        double sum = sum_v[0] + sum_v[1] + sum_v[2] + sum_v[3];
        double sqr_sum = sqr_sum_v[0] + sqr_sum_v[1] + sqr_sum_v[2] + sqr_sum_v[3];
        double mean = sum/nx;
        double norm_sqr_sum = sqr_sum -(2*mean*sum) + (nx*mean*mean);
        double norm_sqrt_sum = sqrt(norm_sqr_sum);

        double4_t norm_sqrt_sum_v = {norm_sqrt_sum, norm_sqrt_sum, norm_sqrt_sum, norm_sqrt_sum};
        double4_t mean_v = {mean, mean, mean, mean};
        double4_t padded_mean_v = {mean, mean, mean, mean};
        if(l == 0) l = 4;
        for (; l < 4; l++) {
            padded_mean_v[l] = 0;
        }
        int s = 0;
        for(; s < row_v_partition_count - 1; s++) {
            vectorized_row[s] = (vectorized_row[s]-mean_v)/norm_sqrt_sum_v;
        }
        vectorized_row[s] = (vectorized_row[s]-padded_mean_v)/norm_sqrt_sum_v;

        vectorized_m.insert(vectorized_m.end(), vectorized_row.begin(), vectorized_row.end());
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double4_t sum = {0, 0, 0, 0};
            for (int k = 0; k < row_v_partition_count; k++) {
                double4_t prod = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
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

    int row_v_partition_count = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m;
    for (int j = 0; j < ny; j++) {

        int n = 0;
        int l = 0;
        double4_t tmp_sum_v = double4_t {0,0,0,0};
        double4_t sum_v = double4_t {0,0,0,0};
        std::vector<double4_t> vectorized_row = std::vector<double4_t> (row_v_partition_count, double4_t {0, 0, 0, 0});
        for (int i = 0; i < nx; i++) {
            double new_dat_p = data[i + nx * j];
            vectorized_row[n][l] = new_dat_p;
            tmp_sum_v[l] = new_dat_p;
            l++;
            if (l>=4) {
                sum_v = sum_v + tmp_sum_v;
                tmp_sum_v = double4_t{0,0,0,0};
                l = 0;
                n++;
            }
        }
        if(l != 0) sum_v = sum_v + tmp_sum_v;

        double row_mean = sum_v[0] + sum_v[1] + sum_v[2] + sum_v[3];
        row_mean = row_mean/nx;
        double4_t vectorized_row_mean = {row_mean, row_mean, row_mean, row_mean};
        double4_t padded_vectorized_row_mean = {row_mean, row_mean, row_mean, row_mean};
        if(l == 0) l = 4;
        for (; l < 4; l++) {
            padded_vectorized_row_mean[l] = 0;
        }
        double4_t sqr_sum_v = double4_t {0,0,0,0};
        int k = 0;
        for (; k < row_v_partition_count - 1; k++) {
            double4_t tmp_value = vectorized_row[k] - vectorized_row_mean;
            vectorized_row[k] = tmp_value;
            sqr_sum_v = sqr_sum_v + tmp_value * tmp_value;

        }
        double4_t tmp_value = vectorized_row[k] - padded_vectorized_row_mean;
        vectorized_row[k] = tmp_value;
        sqr_sum_v = sqr_sum_v + tmp_value * tmp_value;
        double sqrt_sum = sqr_sum_v[0] + sqr_sum_v[1] + sqr_sum_v[2] + sqr_sum_v[3];
        sqrt_sum = sqrt(sqrt_sum);

        for(int s = 0; s < row_v_partition_count; s++) {
            vectorized_row[s] = vectorized_row[s]/sqrt_sum;
        }

        vectorized_m.insert(vectorized_m.end(), vectorized_row.begin(), vectorized_row.end());
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double4_t sum = {0, 0, 0, 0};
            for (int k = 0; k < row_v_partition_count; k++) {
                double4_t prod = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
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