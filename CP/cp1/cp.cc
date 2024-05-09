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

double sum(int start_indx, int end_indx, const float *data) {
    if (start_indx + 1 == end_indx) {
        return double(data[start_indx]) + double(data[end_indx]);
    }
    int mid = start_indx + (end_indx - start_indx)/2;
    if (mid == start_indx) return double(data[mid]);
    return sum(start_indx, mid - 1, data) + sum(mid, end_indx, data);
}

void correlate(int ny, int nx, const float *data, float *result) {

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        double mean = sum(row_indx, row_indx + nx - 1, data);

        mean = mean/nx;
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = data[i + row_indx] - mean;
        }

        double sqr_sum = 0;
        for (int i = 0; i < nx; i++) {
                double val = normalized_data[i + row_indx];
                sqr_sum = sqr_sum + val * val;
        }
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = normalized_data[i + row_indx]/sqrt(sqr_sum);
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            result[i + ny * j] = sum;
        }
    }
}

void correlate(int ny, int nx, const float *data, float *result) {

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        double mean = 0;
        for (int i = 0; i < nx; i++) {
            mean = mean + data[i + row_indx];
        }
        mean = mean/nx;
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = data[i + row_indx] - mean;
        }
        double sqr_sum = 0;
        for (int i = 0; i < nx; i++) {
            sqr_sum = sqr_sum + normalized_data[i + row_indx] * normalized_data[i + row_indx];
        }
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = normalized_data[i + row_indx]/sqrt(sqr_sum);
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            result[i + ny * j] = sum;
        }
    }

}

void correlate(int ny, int nx, const float *data, float *result) {

    int p = 1;

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {
        int row_indx = nx * j;
        double mean_1 = 0;
        double mean_2 = 0;
        double mean_3 = 0;
        double mean_4 = 0;
        for (int i = 0; i < nx/4; i++) {
            mean_1 = mean_1 + data[i + row_indx];
            mean_2 = mean_2 + data[i + nx/4 + row_indx];
            mean_3 = mean_3 + data[i + nx/2 + row_indx];
            mean_4 = mean_4 + data[i + 3*nx/4 + row_indx];
        }

        double mean = (mean_1 + mean_2 + mean_3 + mean_4)/nx;
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = data[i + row_indx] - mean;
        }

        double sqr_sum_1 = 0;
        double sqr_sum_2 = 0;
        double sqr_sum_3 = 0;
        double sqr_sum_4 = 0;
        for (int i = 0; i < nx/4; i++) {
                double val_1 = normalized_data[i + row_indx];
                sqr_sum_1 = sqr_sum_1 + val_1 * val_1;
                double val_2 = normalized_data[i + nx/4 + row_indx];
                sqr_sum_2 = sqr_sum_2 + val_2 * val_2;
                double val_3 = normalized_data[i + nx/2 + row_indx];
                sqr_sum_3 = sqr_sum_3 + val_3 * val_3;
                double val_4 = normalized_data[i + 3*nx/4 + row_indx];
                sqr_sum_4 = sqr_sum_4 + val_4 * val_4;
        }
        double sqr_sum = sqr_sum_1 + sqr_sum_2 + sqr_sum_3 + sqr_sum_4;
        for (int i = 0; i < nx; i++) {
            normalized_data[i + row_indx] = normalized_data[i + row_indx]/sqrt(sqr_sum);
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            result[i + ny * j] = sum;
        }
    }

}

void correlate(int ny, int nx, const float *data, float *result) {

    int p = 4;

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {

        int row_indx = nx * j;
        std::vector<double> sum = std::vector<double> (p, 0);
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                sum[m] = sum[m] + data[m + i*p + row_indx];
            }
        }
        double mean = 0;
        for (int i = 0; i < p; i++) {
            mean = mean + sum[i];
        }
        mean = mean/nx;



        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                normalized_data[m + i*p + row_indx] = data[m + i*p + row_indx] - mean;
            }
        }


        sum = std::vector<double> (p, 0);
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                double val = normalized_data[m + i*p + row_indx];
                sum[m] = sum[m] + val * val;
            }
        }
        double sqr_sum = 0;
        for (int i = 0; i < p; i++) {
            sqr_sum = sqr_sum + sum[i];
        }



        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
            normalized_data[m + i*p + row_indx] = normalized_data[m + i*p + row_indx]/sqrt(sqr_sum);
            }
        }
    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0;
            for (int k = 0; k < nx; k++) {
                sum = sum + normalized_data[k + nx * j] * normalized_data[k + nx * i];
            }
            result[i + ny * j] = sum;
        }
    }
}

void correlate(int ny, int nx, const float *data, float *result) {

    int p = 10;

    double part_1_total = 0;
    double part_2_total = 0;
    double part_3_total = 0;
    double part_4_total = 0;

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {


        auto start = std::chrono::steady_clock::now();
        int row_indx = nx * j;
        std::vector<double> sum = std::vector<double> (p, 0);
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                sum[m] = sum[m] + data[m + i*p + row_indx];
            }
        }
        double mean = 0;
        for (int i = 0; i < p; i++) {
            mean = mean + sum[i];
        }
        mean = mean/nx;
        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        std::chrono::steady_clock::duration diff = end - start;
        part_1_total = part_1_total + double(diff.count()) * std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;



        start = std::chrono::steady_clock::now();
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                normalized_data[m + i*p + row_indx] = data[m + i*p + row_indx] - mean;
            }
        }
        end = std::chrono::steady_clock::now();
        diff = end - start;
        part_2_total = part_2_total + double(diff.count()) * std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;



        start = std::chrono::steady_clock::now();
        sum = std::vector<double> (p, 0);
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                double val = normalized_data[m + i*p + row_indx];
                sum[m] = sum[m] + val * val;
            }
        }
        double sqr_sum = 0;
        for (int i = 0; i < p; i++) {
            sqr_sum = sqr_sum + sum[i];
        }
        end = std::chrono::steady_clock::now();
        diff = end - start;
        part_3_total = part_3_total + double(diff.count()) * std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;



        start = std::chrono::steady_clock::now();
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
            normalized_data[m + i*p + row_indx] = normalized_data[m + i*p + row_indx]/sqrt(sqr_sum);
            }
        }
        end = std::chrono::steady_clock::now();
        diff = end - start;
        part_4_total = part_4_total + double(diff.count()) * std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;


    }



    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            auto sum = std::vector<double> (p, 0);
            for (int k = 0; k < ceil(float(nx)/float(p)); k++) {
                for (int m = 0; m < p && m + k*p < nx; m++) {
                    sum[m] = sum[m] + normalized_data[m + k*p + nx * j] * normalized_data[m + k*p + nx * i];
                }
            }
            for(int m = 0; m < p; m++) {
                result[i + ny * j] = result[i + ny * j] + sum[m];
            }
        }
    }
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::chrono::steady_clock::duration diff = end - start;
    double nseconds = double(diff.count()) * std::chrono::steady_clock::period::num / std::chrono::steady_clock::period::den;
    std::cout << "Part 1 time: " << part_1_total << " ns" << std::endl;
    std::cout << "Part 2 time: " << part_2_total << " ns" << std::endl;
    std::cout << "Part 3 time: " << part_3_total << " ns" << std::endl;
    std::cout << "Part 4 time: " << part_4_total << " ns" << std::endl;
    std::cout << "Part 5 time: " << nseconds << " ns" << std::endl;


}

void correlate(int ny, int nx, const float *data, float *result) {

    int p = 25;

    std::cout << "Partition in use: " << p;

    std::vector<double> normalized_data = std::vector<double> (nx * ny, 0);
    for (int j = 0; j < ny; j++) {

        int row_indx = nx * j;
        std::vector<double> sum = std::vector<double> (p, 0);
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                sum[m] = sum[m] + data[m + i*p + row_indx];
            }
        }
        double mean = 0;
        for (int i = 0; i < p; i++) {
            mean = mean + sum[i];
        }
        mean = mean/nx;

        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                normalized_data[m + i*p + row_indx] = data[m + i*p + row_indx] - mean;
            }
        }

        sum = std::vector<double> (p, 0);
        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
                double val = normalized_data[m + i*p + row_indx];
                sum[m] = sum[m] + val * val;
            }
        }
        double sqr_sum = 0;
        for (int i = 0; i < p; i++) {
            sqr_sum = sqr_sum + sum[i];
        }

        for (int i = 0; i < ceil(float(nx)/float(p)); i++) {
            for (int m = 0; m < p && m + i*p < nx; m++) {
            normalized_data[m + i*p + row_indx] = normalized_data[m + i*p + row_indx]/sqrt(sqr_sum);
            }
        }

    }

    for (int i = 0; i < ny; i++) {
        for (int j = 0; j <= i; j++) {
            std::vector<double> sum = std::vector<double> (p, 0);
            for (int k = 0; k < ceil(float(nx)/float(p)); k++) {
                for (int m = 0; m < p && m + k*p < nx; m++) {
                    sum[m] = sum[m] + normalized_data[m + k*p + nx * j] * normalized_data[m + k*p + nx * i];
                }
            }
            double final_sum = 0;
            for(int m = 0; m < p; m++) {
                final_sum = final_sum + sum[m];
            }
            result[i + ny * j] = final_sum;
        }
    }
}