#include <math.h>
#include <vector>

#include <stdlib.h>
#include <chrono>

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

    auto t1 = std::chrono::high_resolution_clock::now();
    // Ready constants:
    double dividend_by_size = 1/double(nx);

    int row_v_partition_count = ceil(double(nx)/double(4));
    std::vector<double4_t> vectorized_m;
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

        vectorized_m.insert(vectorized_m.end(), vectorized_row.begin(), vectorized_row.end());
    }
    auto t2 = std::chrono::high_resolution_clock::now();


    auto t3 = std::chrono::high_resolution_clock::now();
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
    auto t4 = std::chrono::high_resolution_clock::now();

    auto ns_int_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    //auto ns_int_2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);
    std::cout << ns_int_1.count() << std::endl;
    //std::cout << ns_int_2.count() << std::endl;


}

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