#include <math.h>
#include <vector>

/*
This is the function you need to implement. Quick reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in data[x + y*nx]
- correlation between rows i and row j has to be stored in result[i + j*ny]
- only parts with 0 <= j <= i < ny need to be filled
*/

typedef float float16_t
__attribute__ ((vector_size (16 * sizeof(float))));

typedef double double16_t
__attribute__ ((vector_size (16 * sizeof(double))));

void correlate(int ny, int nx, const float *data, float *result) {

    // Ready constants:
    double dividend_by_size = 1/double(nx);
    int row_v_partition_count = ceil(double(nx)/double(16));

    std::vector<float16_t> vectorized_m = std::vector<float16_t> (ny * row_v_partition_count);
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {

        int n = 0;
        int l = 0;
        double16_t tmp_sum_v = double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        double16_t sum_v = double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        double16_t tmp_sqr_sum_v = double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        double16_t sqr_sum_v = double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        std::vector<double16_t> vectorized_row = std::vector<double16_t> (row_v_partition_count, double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
        for (int i = 0; i < nx; i++) {
            double new_dat_p = data[i + nx * j];
            vectorized_row[n][l] = new_dat_p;
            tmp_sum_v[l] = new_dat_p;
            tmp_sqr_sum_v[l] = new_dat_p * new_dat_p;
            l++;
            if (l>=16) {
                sum_v = sum_v + tmp_sum_v;
                sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
                tmp_sum_v = double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};         // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                tmp_sqr_sum_v = double16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};     // Tämä tosin tarvitsee tehdä vain viimeisessä, joten tässä on vielä parannettavaa
                l = 0;
                n++;
            }
        }
        if(l != 0) {
            sum_v = sum_v + tmp_sum_v;
            sqr_sum_v = sqr_sum_v + tmp_sqr_sum_v;
        }


        double sum = 0;
        double sqr_sum = 0;
        for (int k = 0; k < 16; k++) {
            sum = sum + sum_v[k];
            sqr_sum = sqr_sum + sqr_sum_v[k];
        }

        double mean = sum * dividend_by_size;

        double second_term = -2 * mean;
        double third_term = mean*mean;
        second_term = second_term * sum;
        third_term = third_term * nx;
        sqr_sum = sqr_sum + second_term;
        sqr_sum = sqr_sum + third_term;
        sqr_sum = sqrt(sqr_sum);
        sqr_sum  = 1/sqr_sum;

        double16_t norm_sqrt_sum_v = {sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum, sqr_sum};
        double16_t mean_v = {mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean};
        double16_t padded_mean_v = {mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean, mean};
        if(l == 0) l = 16;
        for (; l < 16; l++) {
            padded_mean_v[l] = 0;
        }
        int s = 0;
        for(; s < row_v_partition_count - 1; s++) {
            double16_t temp = vectorized_row[s]-mean_v;
            vectorized_row[s] = temp * norm_sqrt_sum_v;
        }
        double16_t temp = vectorized_row[s]-padded_mean_v;
        vectorized_row[s] = temp * norm_sqrt_sum_v;

        for (int t = 0; t<row_v_partition_count; t++) {
            for (int k = 0; k<16; k++) {
                vectorized_m[t + j*row_v_partition_count][k] = float(vectorized_row[t][k]);
            }
        }

    }


    int block_size = 4;
    int new_partition_c = ceil(double(row_v_partition_count)/double(4));
    // Padding:
    int padding_size = 0;
    if ((ny % block_size) != 0) padding_size = block_size - (ny % block_size);
    std::vector<float16_t> padding = std::vector<float16_t> (padding_size * row_v_partition_count, float16_t {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
    vectorized_m.insert(vectorized_m.end(), padding.begin(), padding.end());
    int b_ny = ny + padding_size;
    int b_lim = b_ny/block_size;

    #pragma omp parallel for schedule(static,1)
    for (int ib = 0; ib < b_lim; ib++) {
        int i = block_size * ib;
        for (int jb = 0; jb < b_lim; jb++) {
            int j = block_size * jb;

            if (j > i) break;
            float16_t sum_v_00 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_10 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_20 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_30 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_01 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_11 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_21 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_31 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_02 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_12 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_22 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_32 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_03 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_13 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_23 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            float16_t sum_v_33 = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            int k = 0;
            for (; k < new_partition_c - 1; k++) {
                float16_t prod_00_0 = vectorized_m[4 * k + j * row_v_partition_count] * vectorized_m[4 * k + i * row_v_partition_count];
                float16_t prod_00_1 = vectorized_m[1 + 4 * k + j * row_v_partition_count] * vectorized_m[1 + 4 * k + i * row_v_partition_count];
                float16_t prod_00_2 = vectorized_m[2 + 4 * k + j * row_v_partition_count] * vectorized_m[2 + 4 * k + i * row_v_partition_count];
                float16_t prod_00_3 = vectorized_m[3 + 4 * k + j * row_v_partition_count] * vectorized_m[3 + 4 * k + i * row_v_partition_count];

                float16_t prod_10_0 = vectorized_m[4 * k + j * row_v_partition_count] * vectorized_m[4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_10_1 = vectorized_m[1 + 4 * k + j * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_10_2 = vectorized_m[2 + 4 * k + j * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_10_3 = vectorized_m[3 + 4 * k + j * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 1) * row_v_partition_count];

                float16_t prod_20_0 = vectorized_m[4 * k + j * row_v_partition_count] * vectorized_m[4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_20_1 = vectorized_m[1 + 4 * k + j * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_20_2 = vectorized_m[2 + 4 * k + j * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_20_3 = vectorized_m[3 + 4 * k + j * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 2)* row_v_partition_count];

                float16_t prod_30_0 = vectorized_m[4 * k + j * row_v_partition_count] * vectorized_m[4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_30_1 = vectorized_m[1 + 4 * k + j * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_30_2 = vectorized_m[2 + 4 * k + j * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_30_3 = vectorized_m[3 + 4 * k + j * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 3) * row_v_partition_count];

                float16_t prod_01_0 = vectorized_m[4 * k + (j + 1) * row_v_partition_count] * vectorized_m[4 * k + i * row_v_partition_count];
                float16_t prod_01_1 = vectorized_m[1 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[1 + 4 * k + i * row_v_partition_count];
                float16_t prod_01_2 = vectorized_m[2 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[2 + 4 * k + i * row_v_partition_count];
                float16_t prod_01_3 = vectorized_m[3 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[3 + 4 * k + i * row_v_partition_count];

                float16_t prod_11_0 = vectorized_m[4 * k + (j + 1) * row_v_partition_count] * vectorized_m[4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_11_1 = vectorized_m[1 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_11_2 = vectorized_m[2 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_11_3 = vectorized_m[3 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 1) * row_v_partition_count];

                float16_t prod_21_0 = vectorized_m[4 * k + (j + 1) * row_v_partition_count] * vectorized_m[4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_21_1 = vectorized_m[1 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_21_2 = vectorized_m[2 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_21_3 = vectorized_m[3 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 2)* row_v_partition_count];

                float16_t prod_31_0 = vectorized_m[4 * k + (j + 1) * row_v_partition_count] * vectorized_m[4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_31_1 = vectorized_m[1 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_31_2 = vectorized_m[2 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_31_3 = vectorized_m[3 + 4 * k + (j + 1) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 3) * row_v_partition_count];

                float16_t prod_02_0 = vectorized_m[4 * k + (j + 2) * row_v_partition_count] * vectorized_m[4 * k + i * row_v_partition_count];
                float16_t prod_02_1 = vectorized_m[1 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[1 + 4 * k + i * row_v_partition_count];
                float16_t prod_02_2 = vectorized_m[2 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[2 + 4 * k + i * row_v_partition_count];
                float16_t prod_02_3 = vectorized_m[3 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[3 + 4 * k + i * row_v_partition_count];

                float16_t prod_12_0 = vectorized_m[4 * k + (j + 2) * row_v_partition_count] * vectorized_m[4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_12_1 = vectorized_m[1 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_12_2 = vectorized_m[2 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_12_3 = vectorized_m[3 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 1) * row_v_partition_count];

                float16_t prod_22_0 = vectorized_m[4 * k + (j + 2) * row_v_partition_count] * vectorized_m[4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_22_1 = vectorized_m[1 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_22_2 = vectorized_m[2 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_22_3 = vectorized_m[3 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 2)* row_v_partition_count];

                float16_t prod_32_0 = vectorized_m[4 * k + (j + 2) * row_v_partition_count] * vectorized_m[4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_32_1 = vectorized_m[1 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_32_2 = vectorized_m[2 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_32_3 = vectorized_m[3 + 4 * k + (j + 2) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 3) * row_v_partition_count];

                float16_t prod_03_0 = vectorized_m[4 * k + (j + 3) * row_v_partition_count] * vectorized_m[4 * k + i * row_v_partition_count];
                float16_t prod_03_1 = vectorized_m[1 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[1 + 4 * k + i * row_v_partition_count];
                float16_t prod_03_2 = vectorized_m[2 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[2 + 4 * k + i * row_v_partition_count];
                float16_t prod_03_3 = vectorized_m[3 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[3 + 4 * k + i * row_v_partition_count];

                float16_t prod_13_0 = vectorized_m[4 * k + (j + 3) * row_v_partition_count] * vectorized_m[4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_13_1 = vectorized_m[1 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_13_2 = vectorized_m[2 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 1) * row_v_partition_count];
                float16_t prod_13_3 = vectorized_m[3 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 1) * row_v_partition_count];

                float16_t prod_23_0 = vectorized_m[4 * k + (j + 3) * row_v_partition_count] * vectorized_m[4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_23_1 = vectorized_m[1 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_23_2 = vectorized_m[2 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 2)* row_v_partition_count];
                float16_t prod_23_3 = vectorized_m[3 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 2)* row_v_partition_count];

                float16_t prod_33_0 = vectorized_m[4 * k + (j + 3) * row_v_partition_count] * vectorized_m[4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_33_1 = vectorized_m[1 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[1 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_33_2 = vectorized_m[2 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[2 + 4 * k + (i + 3) * row_v_partition_count];
                float16_t prod_33_3 = vectorized_m[3 + 4 * k + (j + 3) * row_v_partition_count] * vectorized_m[3 + 4 * k + (i + 3) * row_v_partition_count];

                sum_v_00 = sum_v_00 + prod_00_0 + prod_00_1 + prod_00_2 + prod_00_3;
                sum_v_10 = sum_v_10 + prod_10_0 + prod_10_1 + prod_10_2 + prod_10_3;
                sum_v_20 = sum_v_20 + prod_20_0 + prod_20_1 + prod_20_2 + prod_20_3;
                sum_v_30 = sum_v_30 + prod_30_0 + prod_30_1 + prod_30_2 + prod_30_3;
                sum_v_01 = sum_v_01 + prod_01_0 + prod_01_1 + prod_01_2 + prod_01_3;
                sum_v_11 = sum_v_11 + prod_11_0 + prod_11_1 + prod_11_2 + prod_11_3;
                sum_v_21 = sum_v_21 + prod_21_0 + prod_21_1 + prod_21_2 + prod_21_3;
                sum_v_31 = sum_v_31 + prod_31_0 + prod_31_1 + prod_31_2 + prod_31_3;
                sum_v_02 = sum_v_02 + prod_02_0 + prod_02_1 + prod_02_2 + prod_02_3;
                sum_v_12 = sum_v_12 + prod_12_0 + prod_12_1 + prod_12_2 + prod_12_3;
                sum_v_22 = sum_v_22 + prod_22_0 + prod_22_1 + prod_22_2 + prod_22_3;
                sum_v_32 = sum_v_32 + prod_32_0 + prod_32_1 + prod_32_2 + prod_32_3;
                sum_v_03 = sum_v_03 + prod_03_0 + prod_03_1 + prod_03_2 + prod_03_3;
                sum_v_13 = sum_v_13 + prod_13_0 + prod_13_1 + prod_13_2 + prod_13_3;
                sum_v_23 = sum_v_23 + prod_23_0 + prod_23_1 + prod_23_2 + prod_23_3;
                sum_v_33 = sum_v_33 + prod_33_0 + prod_33_1 + prod_33_2 + prod_33_3;

            }
            k = k * 4;
            for (; k < row_v_partition_count; k++) {
                float16_t prod_00 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                float16_t prod_10 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                float16_t prod_20 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                float16_t prod_30 = vectorized_m[k + j * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                float16_t prod_01 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                float16_t prod_11 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                float16_t prod_21 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                float16_t prod_31 = vectorized_m[k + (j + 1) * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                float16_t prod_02 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                float16_t prod_12 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                float16_t prod_22 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                float16_t prod_32 = vectorized_m[k + (j + 2) * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
                float16_t prod_03 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + i * row_v_partition_count];
                float16_t prod_13 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + (i + 1) * row_v_partition_count];
                float16_t prod_23 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + (i + 2)* row_v_partition_count];
                float16_t prod_33 = vectorized_m[k + (j + 3) * row_v_partition_count] * vectorized_m[k + (i + 3) * row_v_partition_count];
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

            float sum_00 = 0;
            float sum_10 = 0;
            float sum_20 = 0;
            float sum_30 = 0;
            float sum_01 = 0;
            float sum_11 = 0;
            float sum_21 = 0;
            float sum_31 = 0;
            float sum_02 = 0;
            float sum_12 = 0;
            float sum_22 = 0;
            float sum_32 = 0;
            float sum_03 = 0;
            float sum_13 = 0;
            float sum_23 = 0;
            float sum_33 = 0;
            for (int k = 0; k < 16; k++) {
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

            if (ib < b_lim - 1 && jb < b_lim - 1) {
                result[i + ny * j] = sum_00;
                result[(i + 1) + ny * j] = sum_10;
                result[(i + 2) + ny * j] = sum_20;
                result[(i + 3) + ny * j] = sum_30;
                result[i + ny * (j + 1)] = sum_01;
                result[(i + 1) + ny * (j + 1)] = sum_11;
                result[(i + 2) + ny * (j + 1)] = sum_21;
                result[(i + 3) + ny * (j + 1)] = sum_31;
                result[i + ny * (j + 2)] = sum_02;
                result[(i + 1) + ny * (j + 2)] = sum_12;
                result[(i + 2) + ny * (j + 2)] = sum_22;
                result[(i + 3) + ny * (j + 2)] = sum_32;
                result[i + ny * (j + 3)] = sum_03;
                result[(i + 1) + ny * (j + 3)] = sum_13;
                result[(i + 2) + ny * (j + 3)] = sum_23;
                result[(i + 3) + ny * (j + 3)] = sum_33;
            } else {
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
            }

        }

    }

}