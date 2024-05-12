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

    double delay_0 = 0;
    std::vector<double8_t> vectorized_m = std::vector<double8_t> (ny * row_v_partition_count);
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {

        auto t0 = std::chrono::high_resolution_clock::now();
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
        auto t00 = std::chrono::high_resolution_clock::now();
        auto diff_0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t00 - t0);
        delay_0 = delay_0 + diff_0.count();
    }



    int block_size = 4;
    // Padding:
    int padding_size = 0;
    if ((ny % block_size) != 0) padding_size = block_size - (ny % block_size);
    std::vector<double8_t> padding = std::vector<double8_t> (padding_size * row_v_partition_count, double8_t {0,0,0,0,0,0,0,0});
    vectorized_m.insert(vectorized_m.end(), padding.begin(), padding.end());
    int b_ny = ny + padding_size;
    int b_lim = b_ny/block_size;

    std::cout << "nx: " << nx << ", ny: " << ny << std::endl;
    std::cout << "Padding size: " << padding_size << std::endl;
    std::cout << "b_ny: " << b_ny << std::endl;
    std::cout << "b_lim: " << b_lim << std::endl;

    double delay_1 = 0;
    double delay_2 = 0;
    double delay_3 = 0;
    #pragma omp parallel for schedule(static,1)
    for (int ib = 0; ib < b_lim; ib++) {
        int i = block_size * ib;
        for (int jb = 0; jb < b_lim; jb++) {
            int j = block_size * jb;

            auto t1 = std::chrono::high_resolution_clock::now();
            if (j > i) break;
            double8_t sum_v_00 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_10 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_20 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_30 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_01 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_11 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_21 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_31 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_02 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_12 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_22 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_32 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_03 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_13 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_23 = {0, 0, 0, 0, 0, 0, 0, 0};
            double8_t sum_v_33 = {0, 0, 0, 0, 0, 0, 0, 0};
            for (int k = 0; k < row_v_partition_count; k++) {
                double8_t prod_00 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                double8_t prod_10 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                double8_t prod_20 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                double8_t prod_30 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                double8_t prod_01 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                double8_t prod_11 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                double8_t prod_21 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                double8_t prod_31 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                double8_t prod_02 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                double8_t prod_12 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                double8_t prod_22 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                double8_t prod_32 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                double8_t prod_03 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                double8_t prod_13 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                double8_t prod_23 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                double8_t prod_33 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                sum_v_00 = sum_v_00 + prod_00;
                sum_v_10 = sum_v_10 + prod_10;
                sum_v_20 = sum_v_20 + prod_20;
                sum_v_30 = sum_v_30 + prod_30;
                sum_v_01 = sum_v_01 + prod_01;
                sum_v_11 = sum_v_11 + prod_11;
                sum_v_21 = sum_v_21 + prod_21;
                sum_v_31 = sum_v_31 + prod_31;
                sum_v_02 = sum_v_02 + prod_02;
                sum_v_12 = sum_v_12 + prod_12;
                sum_v_22 = sum_v_22 + prod_22;
                sum_v_32 = sum_v_32 + prod_32;
                sum_v_03 = sum_v_03 + prod_03;
                sum_v_13 = sum_v_13 + prod_13;
                sum_v_23 = sum_v_23 + prod_23;
                sum_v_33 = sum_v_33 + prod_33;
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            auto diff_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
            delay_1 = delay_1 + diff_1.count();

            auto t3 = std::chrono::high_resolution_clock::now();
            double sum_00 = 0;
            double sum_10 = 0;
            double sum_20 = 0;
            double sum_30 = 0;
            double sum_01 = 0;
            double sum_11 = 0;
            double sum_21 = 0;
            double sum_31 = 0;
            double sum_02 = 0;
            double sum_12 = 0;
            double sum_22 = 0;
            double sum_32 = 0;
            double sum_03 = 0;
            double sum_13 = 0;
            double sum_23 = 0;
            double sum_33 = 0;
            for (int k = 0; k < 8; k++) {
                sum_00 = sum_00 + sum_v_00[k];
                sum_10 = sum_10 + sum_v_10[k];
                sum_20 = sum_20 + sum_v_20[k];
                sum_30 = sum_30 + sum_v_30[k];
                sum_01 = sum_01 + sum_v_01[k];
                sum_11 = sum_11 + sum_v_11[k];
                sum_21 = sum_21 + sum_v_21[k];
                sum_31 = sum_31 + sum_v_31[k];
                sum_02 = sum_02 + sum_v_02[k];
                sum_12 = sum_12 + sum_v_12[k];
                sum_22 = sum_22 + sum_v_22[k];
                sum_32 = sum_32 + sum_v_32[k];
                sum_03 = sum_03 + sum_v_03[k];
                sum_13 = sum_13 + sum_v_13[k];
                sum_23 = sum_23 + sum_v_23[k];
                sum_33 = sum_33 + sum_v_33[k];
            }
            auto t4 = std::chrono::high_resolution_clock::now();
            auto diff_2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);
            delay_2 = delay_2 + diff_2.count();

            auto t5 = std::chrono::high_resolution_clock::now();
            result[i + ny * j] = sum_00;
            if (i + 1 < ny) result[(i + 1) + ny * j] = sum_10;
            if (i + 2 < ny) result[(i + 2) + ny * j] = sum_20;
            if (i + 3 < ny) result[(i + 3) + ny * j] = sum_30;
            if (j + 1 < ny) {
                result[i + ny * (j + 1)] = sum_01;
                if (i + 1 < ny) result[(i + 1) + ny * (j + 1)] = sum_11;
                if (i + 2 < ny) result[(i + 2) + ny * (j + 1)] = sum_21;
                if (i + 3 < ny)  result[(i + 3) + ny * (j + 1)] = sum_31;
            }
            if (j + 2 < ny) {
                result[i + ny * (j + 2)] = sum_02;
                if (i + 1 < ny) result[(i + 1) + ny * (j + 2)] = sum_12;
                if (i + 2 < ny) result[(i + 2) + ny * (j + 2)] = sum_22;
                if (i + 3 < ny) result[(i + 3) + ny * (j + 2)] = sum_32;
            }
            if (j + 3 < ny) {
                result[i + ny * (j + 3)] = sum_03;
                if (i + 1 < ny) result[(i + 1) + ny * (j + 3)] = sum_13;
                if (i + 2 < ny) result[(i + 2) + ny * (j + 3)] = sum_23;
                if (i + 3 < ny) result[(i + 3) + ny * (j + 3)] = sum_33;
            }
            auto t6 = std::chrono::high_resolution_clock::now();
            auto diff_3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5);
            delay_3 = delay_3 + diff_3.count();
        }

    }
    std::cout << "Pre processing: " << delay_0 << std::endl;
    std::cout << "First part: " << delay_1 << std::endl;
    std::cout << "Second part: " << delay_2 << std::endl;
    std::cout << "Third part: " << delay_3 << std::endl;

}

/*
void correlate(int ny, int nx, const float *data, float *result) {

    // Delay measurements:
    double delay_0 = 0;
    double delay_1 = 0;
    double delay_2 = 0;
    double delay_3 = 0;

    // Ready constants part 1:
    double dividend_by_size = 1/double(nx);
    int row_v_partition_count = ceil(double(nx)/double(8));

    std::vector<double8_t> vectorized_m = std::vector<double8_t> (ny * row_v_partition_count);
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {

        auto t0 = std::chrono::high_resolution_clock::now();
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
        auto t00 = std::chrono::high_resolution_clock::now();
        auto diff_0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t00 - t0);
        delay_0 = delay_0 + diff_0.count();
    }


    // Ready constants part 2:
    constexpr int block_size = 4;
    // Padding:
    int padding_size = 0;
    if ((ny % block_size) != 0) padding_size = block_size - (ny % block_size);
    std::vector<double8_t> padding = std::vector<double8_t> (padding_size * row_v_partition_count, double8_t {0,0,0,0,0,0,0,0});
    vectorized_m.insert(vectorized_m.end(), padding.begin(), padding.end());
    int b_ny = ny + padding_size;
    int b_lim = b_ny/block_size;

    #pragma omp parallel for
    for (int ib = 0; ib < b_lim; ib++) {
        int i = block_size * ib;
        for (int jb = 0; jb < b_lim; jb++) {
            int j = block_size * jb;
            if (j > i) break;

            auto t1 = std::chrono::high_resolution_clock::now();
            std::vector<double8_t> block_sum_v = std::vector<double8_t> (block_size * block_size, double8_t {0, 0, 0, 0, 0, 0, 0, 0});
            std::vector<double8_t> block_prod_v = std::vector<double8_t> (block_size * block_size, double8_t {0, 0, 0, 0, 0, 0, 0, 0});
            for (int k = 0; k < row_v_partition_count; k++) {
                for (int y = 0; y < block_size; y++) {
                    for (int x = 0; x < block_size; x++) {
                        block_prod_v[x + 4 * y] = vectorized_m[k + (j +  y) * row_v_partition_count] * vectorized_m[k + (i + x) * row_v_partition_count];
                        block_sum_v[x + 4 * y] = block_sum_v[x + 4 * y] + block_prod_v[x + 4 * y];
                    }
                }
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            auto diff_1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
            delay_1 = delay_1 + diff_1.count();

            auto t3 = std::chrono::high_resolution_clock::now();
            std::vector<double> sum = std::vector<double> (block_size * block_size, 0);
            for (int k = 0; k < 8; k++) {
                for (int y = 0; y < block_size; y++) {
                    for (int x = 0; x < block_size; x++) {
                        sum[x + block_size * y] = sum[x + block_size * y] + block_sum_v[x + block_size * y][k];
                    }
                }
            }
            auto t4 = std::chrono::high_resolution_clock::now();
            auto diff_2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);
            delay_2 = delay_2 + diff_2.count();

            auto t5 = std::chrono::high_resolution_clock::now();
            for (int y = 0; y < block_size; y++) {
                for (int x = 0; x < block_size; x++) {
                    if (i + x < ny && j + y < ny) result[(i + x) + ny * (j + y)] = sum[x + block_size * y];
                }
            }
            auto t6 = std::chrono::high_resolution_clock::now();
            auto diff_3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t6 - t5);
            delay_3 = delay_3 + diff_3.count();
        }

    }
    std::cout << "Pre processing: " << delay_0 << std::endl;
    std::cout << "First part: " << delay_1 << std::endl;
    std::cout << "Second part: " << delay_2 << std::endl;
    std::cout << "Third part: " << delay_3 << std::endl;

}
*/