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

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) { 
            normalized_data[i + nx * j] = normalized_data[i + nx * j]/sqrt(sqrd_sums[j]);
        }
    }

    int p = 3;
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            std::vector<double> sum = std::vector<double> (p, 0);
            int k = 0;
            for (; k < floor(nx/p); k++) {
                for (int m = 0; m < p; m++) {
                    double mat = normalized_data[m + k*p + nx * j];
                    double trans = normalized_data[m + k*p + nx * i];
                    double product = mat * trans;
                    sum[m] = sum[m] + product;
                }
            }
            double f_sum = 0;
            for (k = k * p; k < nx; k++) {
                f_sum = f_sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            for(int m = 0; m < p; m++) {
                f_sum = f_sum + sum[m];
            }
            result[i + ny * j] = f_sum;
        }
    }
}