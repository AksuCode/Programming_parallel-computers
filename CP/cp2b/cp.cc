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

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    #pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        
        std::vector<double> row_means = std::vector<double> (ny, 0);
        for (int i = 0; i < nx; i++) {
            row_means[j] = row_means[j] + data[i + nx * j];
        }
        row_means[j] = row_means[j]/nx;
    
        for (int i = 0; i < nx; i++) {
            normalized_data[i + nx * j] = data[i + nx * j] - row_means[j];
        }

        std::vector<double> sqrd_sums = std::vector<double> (ny, 0);
        for (int i = 0; i < nx; i++) {
            sqrd_sums[j] = sqrd_sums[j] + normalized_data[i + nx * j] * normalized_data[i + nx * j];
        }

        for (int i = 0; i < nx; i++) { 
            normalized_data[i + nx * j] = normalized_data[i + nx * j]/sqrt(sqrd_sums[j]);
        }
    }

    #pragma omp parallel for schedule(static,1)
    for (int i = 0; i < ny; i++) {
        std::vector<double> sum_v = std::vector<double> (ny, 0);
        for (int j = 0; j <= i; j++) {
            for (int k = 0; k < nx; k++) {
                sum_v[j] = sum_v[j] + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            result[i + ny * j] = sum_v[j];
        }
    }

}