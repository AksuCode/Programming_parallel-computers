#include <cfloat>
#include <vector>
#include <math.h>
#include <stdlib.h>

#include <chrono>

#include <iostream>

struct Result {
    int y0 = 0;
    int x0 = 0;
    int y1 = 0;
    int x1 = 0;
    float outer[3] = {0};
    float inner[3] = {0};
};


typedef double double4_t
__attribute__ ((vector_size (4 * sizeof(double))));

// Rinnasteisuus ja tehokkaampi suorakulmion koon muuttaminen. Double 4t
Result segment(int ny, int nx, const float *data) {

    // Lets turn data into an array of pixels. Using double4_t for AVX-2.
    std::vector<double4_t> data_pxv(nx * ny);
    double4_t * data_px = &data_pxv[0];
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int c_indx = 3 * i + 3 * nx * j;
            #pragma omp critical
            for (int n = 0; n < 3; n++) {
                data_px[i + nx * j][n] = double(data[n + c_indx]);
            }
        }
    }

    // Sum of all pixel color component values stored in vector (R,G,B).
    double4_t color_sum_v = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            #pragma omp critical
            color_sum_v = color_sum_v + data_px[i + nx * j];
        }
    }

    //double csqr_sum = 0;    // Squared sum for all pixels and colors
    double4_t csqr_sum = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            csqr_sum = csqr_sum + data_px[i + nx * j] * data_px[i + nx * j];
        }
    }

    // Initialize inside regions
    std::vector<double4_t> initial_inside_numerator_v((nx + 1) * (ny + 1), double4_t {0,0,0,0});
    double4_t * initial_inside_numerator = &initial_inside_numerator_v[0];
    #pragma omp parallel for collapse(2)
    for (int w = 1; w <= nx; w++) {
        for (int h = 1; h <= ny; h++) {
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    initial_inside_numerator[w + nx * h] = initial_inside_numerator[w + nx * h] + data_px[i + nx *j];
                }
            }
        }
    }

    // Save results
    std::pair<double, Result> dummy;
    dummy.first = DBL_MAX;
    std::vector<std::pair<double, Result>> results(nx + 1, dummy);

    // Take time
    double t1 = 0;
    double t2 = 0;
    double t3 = 0;

    // Create rectangle segments of varying size.
    //#pragma omp parallel for schedule(static,1)
    for (int w = 1; w <= nx; w++) {

        // Initialize stuff for the threads:
        // Minimum cost:
        double minimum_cost = DBL_MAX;
        struct Result res;

        for (int h = 1; h <= ny; h++) {

            // Denumerators
            double4_t inside_denumerator = {1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w)};
            double4_t outside_denumerator = double4_t {1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w))};

            // Inside region
            std::vector<double4_t> alloc_inside_numerator_v(nx * ny, double4_t {0,0,0,0});
            double4_t * inside_numerator_v = &alloc_inside_numerator_v[0];
            inside_numerator_v[0] = initial_inside_numerator[w + nx * h];

            for (int j = 0; j < h; j++) {
                for (int i = 1; i <= nx - w; i++) {
                    int c_indx = (i - 1) + nx * j;
                    inside_numerator_v[i] = inside_numerator_v[i - 1] - data_px[c_indx] + data_px[c_indx + w];
                }
            }
            for (int j = 1; j < ny - h) {
                for (int i = 0; i < nx - w; i++) {
                    int c_indx = i + nx * (j - 1);
                    int h_diff = nx * h;
                    inside_numerator_v[i + nx * j] = inside_numerator_v[i + nx * (j - 1)] - data_px[i + nx * (j - 1)] + data_px[c_indx + h_diff];
                }
            }

            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                auto t31 = std::chrono::high_resolution_clock::now();
                for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {

                    int low_right_c_x = w + up_left_c_x;   // Lower right corner x pos

                    auto t21 = std::chrono::high_resolution_clock::now();
                    if (up_left_c_y != 0) {
                        for (int i = up_left_c_x; i < low_right_c_x; i++) {
                            int c_indx = i + nx * (up_left_c_y - 1);
                            int h_diff = nx * h;
                            inside_numerator_v[up_left_c_x] = inside_numerator_v[up_left_c_x] - data_px[c_indx] + data_px[c_indx + h_diff];
                        }
                    }
                    auto t22 = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double, std::milli> inst2 = t22-t21;
                    t2 = t2 + inst2.count();

                    double4_t outside_numerator_v = color_sum_v - inside_numerator_v[up_left_c_x];

                    // Calculating cost:
                    double4_t inside_sqr_num = inside_numerator_v[up_left_c_x] * inside_numerator_v[up_left_c_x];
                    double4_t outside_sqr_num = outside_numerator_v * outside_numerator_v;
                    double4_t cost = csqr_sum -(inside_sqr_num*inside_denumerator) -(outside_sqr_num*outside_denumerator);

                    double final_cost = cost[0] + cost[1] + cost[2];

                    if (final_cost < minimum_cost) {

                        res.x0 = up_left_c_x;
                        res.x1 = up_left_c_x + w;
                        res.y0 = up_left_c_y;
                        res.y1 = low_right_c_y;
                        double4_t in = inside_numerator_v[up_left_c_x]*inside_denumerator;
                        double4_t out = outside_numerator_v*outside_denumerator;
                        res.inner[0] = in[0];
                        res.inner[1] = in[1];
                        res.inner[2] = in[2];
                        res.outer[0] = out[0];
                        res.outer[1] = out[1];
                        res.outer[2] = out[2];

                        minimum_cost = final_cost;

                    }
                }
                auto t32 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli> inst3 = t32-t31;
                t3 = t3 + inst3.count();
            }
        }

        results[w] = std::pair(minimum_cost, res);
    }

    std::cout << "First: " << t1 << std::endl;
    std::cout << "Second: " << t2 << std::endl;
    std::cout << "Third: " << t3 << std::endl;

    int res_indx = 0;
    double minimum_cost = DBL_MAX;

    for(int k = 0; k < nx + 1; k++) {
        std::pair<double,Result> tmp = results[k];
        if (tmp.first < minimum_cost) {
            minimum_cost = tmp.first;
            res_indx = k;
        }
    }

    return results[res_indx].second;
}



/*
typedef double double4_t
__attribute__ ((vector_size (4 * sizeof(double))));

// Rinnasteisuus ja tehokkaampi suorakulmion koon muuttaminen. Double 4t
Result segment(int ny, int nx, const float *data) {

    // Lets turn data into an array of pixels. Using double4_t for AVX-2.
    std::vector<double4_t> data_pxv(nx * ny);
    double4_t * data_px = &data_pxv[0];
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int c_indx = 3 * i + 3 * nx * j;
            #pragma omp critical
            for (int n = 0; n < 3; n++) {
                data_px[i + nx * j][n] = double(data[n + c_indx]);
            }
        }
    }

    // Sum of all pixel color component values stored in vector (R,G,B).
    double4_t color_sum_v = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            #pragma omp critical
            color_sum_v = color_sum_v + data_px[i + nx * j];
        }
    }

    //double csqr_sum = 0;    // Squared sum for all pixels and colors
    double4_t csqr_sum = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            csqr_sum = csqr_sum + data_px[i + nx * j] * data_px[i + nx * j];
        }
    }

    // Initialize inside regions
    std::vector<double4_t> initial_inside_numerator_v((nx + 1) * (ny + 1), double4_t {0,0,0,0});
    double4_t * initial_inside_numerator = &initial_inside_numerator_v[0];
    #pragma omp parallel for collapse(2)
    for (int w = 1; w <= nx; w++) {
        for (int h = 1; h <= ny; h++) {
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    initial_inside_numerator[w + nx * h] = initial_inside_numerator[w + nx * h] + data_px[i + nx *j];
                }
            }
        }
    }

    // Save results
    std::pair<double, Result> dummy;
    dummy.first = DBL_MAX;
    std::vector<std::pair<double, Result>> results(nx + 1, dummy);

    // Take time
    double t1 = 0;
    double t2 = 0;
    double t3 = 0;

    // Create rectangle segments of varying size.
    //#pragma omp parallel for schedule(static,1)
    for (int w = 1; w <= nx; w++) {

        // Initialize stuff for the threads:
        // Minimum cost:
        double minimum_cost = DBL_MAX;
        struct Result res;

        for (int h = 1; h <= ny; h++) {

            // Denumerators
            double4_t inside_denumerator = {1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w)};
            double4_t outside_denumerator = double4_t {1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w))};

            // Inside region
            double4_t inside_numerator = initial_inside_numerator[w + nx * h];

            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                auto t11 = std::chrono::high_resolution_clock::now();
                if (up_left_c_y != 0) {
                    // Go through values that are left outside the new region our that are added to the new region after movement in y direction:
                    for (int i = 0; i < w; i++) {
                        int c_indx = i + nx * (up_left_c_y - 1);
                        int h_diff = nx * h;
                        inside_numerator = inside_numerator - data_px[c_indx] + data_px[c_indx + h_diff];
                    }
                }
                auto t12 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli> inst1 = t12-t11;
                t1 = t1 + inst1.count();

                auto t21 = std::chrono::high_resolution_clock::now();
                // Save values for movement in x direction
                std::vector<double4_t> alloc_tmp_inside_numerator_v(nx-w + 1, double4_t {0,0,0,0});
                double4_t * tmp_inside_numerator_v = &alloc_tmp_inside_numerator_v[0];
                for (int j = up_left_c_y; j < low_right_c_y; j++) {
                    for (int up_left_c_x = 1; w + up_left_c_x <= nx; up_left_c_x++) {
                        int c_indx = (up_left_c_x - 1) + nx * j;
                        tmp_inside_numerator_v[up_left_c_x] = tmp_inside_numerator_v[up_left_c_x] - data_px[c_indx] + data_px[c_indx + w];
                    }
                }
                auto t22 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli> inst2 = t22-t21;
                t2 = t2 + inst2.count();

                auto t31 = std::chrono::high_resolution_clock::now();
                double4_t cumulative_in = inside_numerator;
                for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {
                    cumulative_in = cumulative_in + tmp_inside_numerator_v[up_left_c_x];
                    double4_t cumulative_out = color_sum_v - cumulative_in;

                    // Calculating cost:
                    double4_t inside_sqr_num = cumulative_in * cumulative_in;
                    double4_t outside_sqr_num = cumulative_out * cumulative_out;
                    double4_t cost = csqr_sum -(inside_sqr_num*inside_denumerator) -(outside_sqr_num*outside_denumerator);

                    double final_cost = cost[0] + cost[1] + cost[2];

                    if (final_cost < minimum_cost) {

                        res.x0 = up_left_c_x;
                        res.x1 = up_left_c_x + w;
                        res.y0 = up_left_c_y;
                        res.y1 = low_right_c_y;
                        double4_t in = cumulative_in*inside_denumerator;
                        double4_t out = cumulative_out*outside_denumerator;
                        res.inner[0] = in[0];
                        res.inner[1] = in[1];
                        res.inner[2] = in[2];
                        res.outer[0] = out[0];
                        res.outer[1] = out[1];
                        res.outer[2] = out[2];

                        minimum_cost = final_cost;

                    }
                }
                auto t32 = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double, std::milli> inst3 = t32-t31;
                t3 = t3 + inst3.count();
            }
        }

        results[w] = std::pair(minimum_cost, res);
    }

    std::cout << "First: " << t1 << std::endl;
    std::cout << "Second: " << t2 << std::endl;
    std::cout << "Third: " << t3 << std::endl;

    int res_indx = 0;
    double minimum_cost = DBL_MAX;

    for(int k = 0; k < nx + 1; k++) {
        std::pair<double,Result> tmp = results[k];
        if (tmp.first < minimum_cost) {
            minimum_cost = tmp.first;
            res_indx = k;
        }
    }

    return results[res_indx].second;
}
*/

/*
typedef double double8_t
__attribute__ ((vector_size (8 * sizeof(double))));

// Rinnasteisuus ja tehokkaampi suorakulmion koon muuttaminen. Double 4t
Result segment(int ny, int nx, const float *data) {

    // ILP
    constexpr int nb = 4;
    int v_row_num = (nx + 8 - 1)/8;
    int na = (v_row_num + nb - 1) / nb;
    int nab = na*nb;

    // Vectorization
    std::vector<double8_t> data_R_v(nab * ny);
    std::vector<double8_t> data_G_v(nab * ny);
    std::vector<double8_t> data_B_v(nab * ny);
    double8_t * data_R = &data_R_v[0];
    double8_t * data_G = &data_G_v[0];
    double8_t * data_B = &data_B_v[0];
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int c_indx = 3 * i + 3 * nx * j;
            #pragma omp critical
            for (int n = 0; n < 8; n++) {
                data_px[i + nab * j][n] = double(data[n + c_indx]);
            }
        }
    }

    // Sum of all pixel color component values stored in vector (R,G,B).
    double4_t color_sum_v = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            #pragma omp critical
            color_sum_v = color_sum_v + data_px[i + nx * j];
        }
    }

    //double csqr_sum = 0;    // Squared sum for all pixels and colors
    double4_t csqr_sum = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            csqr_sum = csqr_sum + data_px[i + nx * j] * data_px[i + nx * j];
        }
    }


    std::pair<double, Result> dummy;
    dummy.first = DBL_MAX;
    std::vector<std::pair<double, Result>> results(nx + 1, dummy);

    // Create rectangle segments of varying size.
    #pragma omp parallel for schedule(static,1)
    for (int w = 1; w <= nx; w++) {

        // Inside region (w)
        double4_t initial_inside_numerator = {0,0,0,0};  // Inside segment (R, G, B)
        for (int i = 0; i < w; i++) {
            initial_inside_numerator = initial_inside_numerator + data_px[i];
        }


        // Initialize stuff for the threads:
        // Minimum cost:
        double minimum_cost = DBL_MAX;
        struct Result res;

        for (int h = 1; h <= ny; h++) {

            if (h != 1) {
                for (int i = 0; i < w; i++) {
                    initial_inside_numerator = initial_inside_numerator + data_px[i + nx * (h - 1)];
                }
            }

            // Denumerators
            double4_t inside_denumerator = {1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w)};
            double4_t outside_denumerator = double4_t {1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w))};

            // Numerators
            // Inside region
            double4_t inside_numerator = initial_inside_numerator;
            // Outside region
            double4_t outside_numerator = color_sum_v - inside_numerator;

            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                if (up_left_c_y != 0) {
                    // Go through values that are left outside the new region our that are added to the new region after movement in y direction:
                    for (int i = 0; i < w; i++) {
                        int c_indx = i + nx * (up_left_c_y - 1);
                        int h_diff = nx * h;
                        outside_numerator = outside_numerator + data_px[c_indx] - data_px[c_indx + h_diff];
                        inside_numerator = inside_numerator - data_px[c_indx] + data_px[c_indx + h_diff];
                    }
                }

                // Save values for movement in x direction
                std::vector<double4_t> alloc_tmp_outside_numerator_v(nx-w + 1, double4_t {0,0,0,0});
                double4_t * tmp_outside_numerator_v = &alloc_tmp_outside_numerator_v[0];
                std::vector<double4_t> alloc_tmp_inside_numerator_v(nx-w + 1, double4_t {0,0,0,0});
                double4_t * tmp_inside_numerator_v = &alloc_tmp_inside_numerator_v[0];
                for (int j = up_left_c_y; j < low_right_c_y; j++) {
                    for (int up_left_c_x = 1; w + up_left_c_x <= nx; up_left_c_x++) {
                        int c_indx = (up_left_c_x - 1) + nx * j;
                        tmp_inside_numerator_v[up_left_c_x] = tmp_inside_numerator_v[up_left_c_x] - data_px[c_indx] + data_px[c_indx + w];
                        tmp_outside_numerator_v[up_left_c_x] = tmp_outside_numerator_v[up_left_c_x] + data_px[c_indx] - data_px[c_indx + w]; 
                    }
                }

                double4_t cumulative_in = inside_numerator;
                double4_t cumulative_out = outside_numerator;
                for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {
                    cumulative_in = cumulative_in + tmp_inside_numerator_v[up_left_c_x];
                    cumulative_out = cumulative_out + tmp_outside_numerator_v[up_left_c_x];

                    // Calculating cost:
                    double4_t inside_sqr_num = cumulative_in * cumulative_in;
                    double4_t outside_sqr_num = cumulative_out * cumulative_out;
                    double4_t cost = csqr_sum -(inside_sqr_num*inside_denumerator) -(outside_sqr_num*outside_denumerator);

                    double final_cost = cost[0] + cost[1] + cost[2];

                    if (final_cost < minimum_cost) {

                        res.x0 = up_left_c_x;
                        res.x1 = up_left_c_x + w;
                        res.y0 = up_left_c_y;
                        res.y1 = low_right_c_y;
                        double4_t in = cumulative_in*inside_denumerator;
                        double4_t out = cumulative_out*outside_denumerator;
                        res.inner[0] = in[0];
                        res.inner[1] = in[1];
                        res.inner[2] = in[2];
                        res.outer[0] = out[0];
                        res.outer[1] = out[1];
                        res.outer[2] = out[2];

                        minimum_cost = final_cost;

                    }
                }
            }
        }

        results[w] = std::pair(minimum_cost, res);
    }

    int res_indx = 0;
    double minimum_cost = DBL_MAX;

    for(int k = 0; k < nx + 1; k++) {
        std::pair<double,Result> tmp = results[k];
        if (tmp.first < minimum_cost) {
            minimum_cost = tmp.first;
            res_indx = k;
        }
    }

    return results[res_indx].second;
}
*/

/*
typedef double double4_t
__attribute__ ((vector_size (4 * sizeof(double))));

// Rinnasteisuus ja tehokkaampi suorakulmion koon muuttaminen. Double 4t
Result segment(int ny, int nx, const float *data) {

    // Lets turn data into an array of pixels. Using double4_t for AVX-2.
    std::vector<double4_t> data_pxv(nx * ny);
    double4_t * data_px = &data_pxv[0];
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int c_indx = 3 * i + 3 * nx * j;
            #pragma omp critical
            for (int n = 0; n < 3; n++) {
                data_px[i + nx * j][n] = double(data[n + c_indx]);
            }
        }
    }

    // Sum of all pixel color component values stored in vector (R,G,B).
    double4_t color_sum_v = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            #pragma omp critical
            color_sum_v = color_sum_v + data_px[i + nx * j];
        }
    }

    //double csqr_sum = 0;    // Squared sum for all pixels and colors
    double4_t csqr_sum = {0,0,0,0};
    #pragma omp parallel for schedule(static,1)
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            csqr_sum = csqr_sum + data_px[i + nx * j] * data_px[i + nx * j];
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();
    // Initialize inside regions
    std::vector<double4_t> initial_inside_numerator_v((nx + 1) * (ny + 1), double4_t {0,0,0,0});
    double4_t * initial_inside_numerator = &initial_inside_numerator_v[0];
    #pragma omp parallel for collapse(2)
    for (int w = 1; w <= nx; w++) {
        for (int h = 1; h <= ny; h++) {
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    initial_inside_numerator[w + nx * h] = initial_inside_numerator[w + nx * h] + data_px[i + nx *j];
                }
            }
        }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> t2t1 = t2 - t1;
    std::cout << t2t1.count() << std::endl;

    // Save results
    std::pair<double, Result> dummy;
    dummy.first = DBL_MAX;
    std::vector<std::pair<double, Result>> results(nx + 1, dummy);

    auto t3 = std::chrono::high_resolution_clock::now();
    // Create rectangle segments of varying size.
    #pragma omp parallel for schedule(static,1)
    for (int w = 1; w <= nx; w++) {

        // Initialize stuff for the threads:
        // Minimum cost:
        double minimum_cost = DBL_MAX;
        struct Result res;

        for (int h = 1; h <= ny; h++) {

            // Denumerators
            double4_t inside_denumerator = {1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w), 1.0/double(h * w)};
            double4_t outside_denumerator = double4_t {1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w)), 1.0/(double(nx * ny)-double(h * w))};

            // Numerators
            // Inside region
            double4_t inside_numerator = initial_inside_numerator[w + nx * h];
            // Outside region
            double4_t outside_numerator = color_sum_v - inside_numerator;

            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                if (up_left_c_y != 0) {
                    // Go through values that are left outside the new region our that are added to the new region after movement in y direction:
                    for (int i = 0; i < w; i++) {
                        int c_indx = i + nx * (up_left_c_y - 1);
                        int h_diff = nx * h;
                        outside_numerator = outside_numerator + data_px[c_indx] - data_px[c_indx + h_diff];
                        inside_numerator = inside_numerator - data_px[c_indx] + data_px[c_indx + h_diff];
                    }
                }

                // Save values for movement in x direction
                std::vector<double4_t> alloc_tmp_outside_numerator_v(nx-w + 1, double4_t {0,0,0,0});
                double4_t * tmp_outside_numerator_v = &alloc_tmp_outside_numerator_v[0];
                std::vector<double4_t> alloc_tmp_inside_numerator_v(nx-w + 1, double4_t {0,0,0,0});
                double4_t * tmp_inside_numerator_v = &alloc_tmp_inside_numerator_v[0];
                for (int j = up_left_c_y; j < low_right_c_y; j++) {
                    for (int up_left_c_x = 1; w + up_left_c_x <= nx; up_left_c_x++) {
                        int c_indx = (up_left_c_x - 1) + nx * j;
                        tmp_inside_numerator_v[up_left_c_x] = tmp_inside_numerator_v[up_left_c_x] - data_px[c_indx] + data_px[c_indx + w];
                        tmp_outside_numerator_v[up_left_c_x] = tmp_outside_numerator_v[up_left_c_x] + data_px[c_indx] - data_px[c_indx + w]; 
                    }
                }

                double4_t cumulative_in = inside_numerator;
                double4_t cumulative_out = outside_numerator;
                for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {
                    cumulative_in = cumulative_in + tmp_inside_numerator_v[up_left_c_x];
                    cumulative_out = cumulative_out + tmp_outside_numerator_v[up_left_c_x];

                    // Calculating cost:
                    double4_t inside_sqr_num = cumulative_in * cumulative_in;
                    double4_t outside_sqr_num = cumulative_out * cumulative_out;
                    double4_t cost = csqr_sum -(inside_sqr_num*inside_denumerator) -(outside_sqr_num*outside_denumerator);

                    double final_cost = cost[0] + cost[1] + cost[2];

                    if (final_cost < minimum_cost) {

                        res.x0 = up_left_c_x;
                        res.x1 = up_left_c_x + w;
                        res.y0 = up_left_c_y;
                        res.y1 = low_right_c_y;
                        double4_t in = cumulative_in*inside_denumerator;
                        double4_t out = cumulative_out*outside_denumerator;
                        res.inner[0] = in[0];
                        res.inner[1] = in[1];
                        res.inner[2] = in[2];
                        res.outer[0] = out[0];
                        res.outer[1] = out[1];
                        res.outer[2] = out[2];

                        minimum_cost = final_cost;

                    }
                }
            }
        }

        results[w] = std::pair(minimum_cost, res);
    }
    auto t4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> t4t3 = t4 - t3;
    std::cout << t4t3.count() << std::endl;

    int res_indx = 0;
    double minimum_cost = DBL_MAX;

    for(int k = 0; k < nx + 1; k++) {
        std::pair<double,Result> tmp = results[k];
        if (tmp.first < minimum_cost) {
            minimum_cost = tmp.first;
            res_indx = k;
        }
    }

    return results[res_indx].second;
}
*/

/*
typedef double double8_t
__attribute__ ((vector_size (8 * sizeof(double))));

//Double 8t implementaatio (vaiheessa) ja lineaarisempi lukeminen sisäisimmällä loopilla:
Result segment(int ny, int nx, const float *data) {

    // Lets turn data into an array of pixels. Using double4_t for AVX-2.
    int row_v_num = ceil(double(nx)/double(8));

    // Create arrays for each color
    std::vector<double8_t> data_R(row_v_num * ny);
    std::vector<double8_t> data_G(row_v_num * ny);
    std::vector<double8_t> data_B(row_v_num * ny);
    //#pragma omp parallel for
    for (int y = 0; y < ny; y++) {
        int l = 0;
        int s = 0;
        for (int x = 0; x < nx; x++) {
            int indx = 3 * x + 3 * nx * y;
            data_R[l][s] = data[indx];
            data_G[l][s] = data[1 + indx];
            data_B[l][s] = data[2 + indx];
            if (s >= 8) {
                s=0;
                l++;
            }
        }
    }

    // Sum of all pixel color component values
    double8_t sum_R = {0,0,0,0,0,0,0,0};
    double8_t sum_G = {0,0,0,0,0,0,0,0};
    double8_t sum_B = {0,0,0,0,0,0,0,0};
    // Squared sum for all pixels and colors
    double8_t sqr_sum = {0,0,0,0,0,0,0,0};
    //#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < row_v_num; i++) {
            double8_t sqr_sum_R = data_R[i + row_v_num * j] * data_R[i + row_v_num * j];
            double8_t sqr_sum_G = data_R[i + row_v_num * j] * data_G[i + row_v_num * j];
            double8_t sqr_sum_B = data_R[i + row_v_num * j] * data_B[i + row_v_num * j];
            sum_R = sum_R + data_R[i + row_v_num * j];
            sum_G = sum_G + data_G[i + row_v_num * j];
            sum_B = sum_B + data_B[i + row_v_num * j];
            sqr_sum = sqr_sum + sqr_sum_R;
            sqr_sum = sqr_sum + sqr_sum_G;
            sqr_sum = sqr_sum + sqr_sum_B;
        }
    }

    struct Result res;
    double minimum_cost = DBL_MAX;
    
    // Create rectangle segments of varying size.
    for (int w = 1; w <= nx; w++) {
        for (int h = 1; h <= ny; h++) {

            // Inside region
            double8_t inside_num_R = {0,0,0,0,0,0,0,0};
            double8_t inside_num_G = {0,0,0,0,0,0,0,0};
            double8_t inside_num_B = {0,0,0,0,0,0,0,0};
            int vectorized_w = w/8;
            for (int j = 0; j < h; j++) {
                int i = 0;
                for (; i < vectorized_w; i++) {
                    inside_num_R = inside_num_R + data_R[i + row_v_num * j];
                    inside_num_G = inside_num_G + data_G[i + row_v_num * j];
                    inside_num_B = inside_num_B + data_B[i + row_v_num * j];
                }
                int lim = w - i * 8;
                for (int n = 0; n < lim; n++) {
                    inside_num_R[n] = inside_num_R[n] + data_R[i + row_v_num * j][n];
                    inside_num_G[n] = inside_num_G[n] + data_G[i + row_v_num * j][n];
                    inside_num_B[n] = inside_num_B[n] + data_B[i + row_v_num * j][n];
                }
            }

            // Outside region
            double8_t outside_num_R = sum_R - inside_num_R;
            double8_t outside_num_G = sum_G - inside_num_G;
            double8_t outside_num_B = sum_B - inside_num_B;

            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {
                    
                int low_right_c_x = w + up_left_c_x;

                if (up_left_c_x != 0) {
                    // Go through values that are left outside the new region our that are added to the new region after movement in x direction:
                    for (int j = 0; j < h; j++) {
                        int c_indx = (up_left_c_x - 1) + row_v_num * j;
                        int w_diff = w;
                        for (int n = 0; n < 3; n++) {
                            outside_numerator[n] = outside_numerator[n] + data[n+ c_indx] - data[n + c_indx + w_diff];
                            inside_numerator[n] = inside_numerator[n] - data[n+ c_indx] + data[n + c_indx + w_diff];
                        }
                    } 
                }

                // Save values for movement in y direction
                std::vector<double> tmp_outside_numerator = outside_numerator;
                std::vector<double> tmp_inside_numerator = inside_numerator;

                for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                    int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                    if (up_left_c_y != 0) {
                        // Go through values that are left outside the new region our that are added to the new region after movement in y direction:
                        for (int i = up_left_c_x; i < low_right_c_x; i++) {
                            int c_indx = 3 * i + 3 * nx * (up_left_c_y - 1);
                            int h_diff = 3 * nx * h;
                            for (int n = 0; n < 3; n++) {
                                tmp_outside_numerator[n] = tmp_outside_numerator[n] + data[n + c_indx] - data[n + c_indx + h_diff];
                                tmp_inside_numerator[n] = tmp_inside_numerator[n] - data[n + c_indx] + data[n + c_indx + h_diff];
                            }
                        }
                    }

                    // Calculating cost:
                    int inside_denumerator = h * w;
                    int outside_denumerator = nx * ny - inside_denumerator;
                    double inside_sqr_num = 0;
                    double outside_sqr_num = 0;
                    for (int n = 0; n < 3; n++) {
                        inside_sqr_num = inside_sqr_num + double(tmp_inside_numerator[n]) * double(tmp_inside_numerator[n]);
                        outside_sqr_num = outside_sqr_num + double(tmp_outside_numerator[n]) * double(tmp_outside_numerator[n]);
                    }
                    double cost = csqr_sum -(inside_sqr_num/inside_denumerator) -(outside_sqr_num/outside_denumerator);

                    // Check if a new minimum cost is found:
                    if (cost < minimum_cost) {
                        minimum_cost = cost;
                        res.x0 = up_left_c_x;
                        res.x1 = low_right_c_x;
                        res.y0 = up_left_c_y;
                        res.y1 = low_right_c_y;
                        for (int n = 0; n < 3; n++) {
                            res.inner[n] = tmp_inside_numerator[n]/inside_denumerator;
                            res.outer[n] = tmp_outside_numerator[n]/outside_denumerator;
                        }
                    }
                }
            }
        }
    }

    return res;
}
*/

/*
// Hieman rikki oleva double8t implementatio. Rinnakkainen
Result segment(int ny, int nx, const float *data) {

    // Lets turn data into an array of pixels. Using double4_t for AVX-2.
    int row_v_num = ceil(double(nx)/double(8));

    // Create arrays for each color
    std::vector<double8_t> data_R(row_v_num * ny);
    std::vector<double8_t> data_G(row_v_num * ny);
    std::vector<double8_t> data_B(row_v_num * ny);
    //#pragma omp parallel for
    for (int y = 0; y < ny; y++) {
        int l = 0;
        int s = 0;
        for (int x = 0; x < nx; x++) {
            int indx = 3 * x + 3 * nx * y;
            data_R[l][s] = data[indx];
            data_G[l][s] = data[1 + indx];
            data_B[l][s] = data[2 + indx];
            if (s >= 8) {
                s=0;
                l++;
            }
        }
    }

    // Sum of all pixel color component values
    double8_t sum_R = {0,0,0,0,0,0,0,0};
    double8_t sum_G = {0,0,0,0,0,0,0,0};
    double8_t sum_B = {0,0,0,0,0,0,0,0};
    //#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < row_v_num; i++) {
            sum_R = sum_R + data_R[i + row_v_num * j];
            sum_G = sum_G + data_G[i + row_v_num * j];
            sum_B = sum_B + data_B[i + row_v_num * j];
        }
    }

    // Squared sum for all pixels and colors
    double8_t sqr_sum = {0,0,0,0,0,0,0,0};
    //#pragma omp parallel for
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < row_v_num; i++) {
            double8_t sqr_sum_R = data_R[i + row_v_num * j] * data_R[i + row_v_num * j];
            double8_t sqr_sum_G = data_R[i + row_v_num * j] * data_G[i + row_v_num * j];
            double8_t sqr_sum_B = data_R[i + row_v_num * j] * data_B[i + row_v_num * j];
            sqr_sum = sqr_sum + sqr_sum_R;
            sqr_sum = sqr_sum + sqr_sum_G;
            sqr_sum = sqr_sum + sqr_sum_B;
        }
    }

    std::pair<double, Result> dummy;
    dummy.first = DBL_MAX;
    std::vector<std::pair<double, Result>> results(nx + 1, dummy);
    
    // Create rectangle segments of varying size.
    //#pragma omp parallel for schedule(static,1)
    for (int w = 1; w <= nx; w++) {
        std::vector<std::pair<double, Result>> tmp_results(ny + 1, dummy);
        for (int h = 1; h <= ny; h++) {

            // Inside region
            double8_t inside_num_R = {0,0,0,0,0,0,0,0};
            double8_t inside_num_G = {0,0,0,0,0,0,0,0};
            double8_t inside_num_B = {0,0,0,0,0,0,0,0};
            int vectorized_w = w/8;
            for (int j = 0; j < h; j++) {
                int i = 0;
                for (; i < vectorized_w; i++) {
                    inside_num_R = inside_num_R + data_R[i + row_v_num * j];
                    inside_num_G = inside_num_G + data_G[i + row_v_num * j];
                    inside_num_B = inside_num_B + data_B[i + row_v_num * j];
                }
                int lim = w - i * 8;
                for (int n = 0; n < lim; n++) {
                    inside_num_R[n] = inside_num_R[n] + data_R[i + row_v_num * j][n];
                    inside_num_G[n] = inside_num_G[n] + data_G[i + row_v_num * j][n];
                    inside_num_B[n] = inside_num_B[n] + data_B[i + row_v_num * j][n];
                }
            }

            // Outside region
            double8_t outside_num_R = sum_R - inside_num_R;
            double8_t outside_num_G = sum_G - inside_num_G;
            double8_t outside_num_B = sum_B - inside_num_B;

            int tmp_results_2_it = 0;
            std::vector<std::pair<double, Result>> tmp_results_2((1 + nx) * (1 + ny), dummy);
            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                if (up_left_c_y != 0) {
                    // Go through values that are left outside the new region our that are added to the new region after movement in y direction:
                    int i = 0;
                    for (; i < vectorized_w; i++) {
                        int c_indx = i + row_v_num * (up_left_c_y - 1);
                        int h_diff = row_v_num * h;
                        outside_num_R = outside_num_R + data_R[c_indx] - data_R[c_indx + h_diff];
                        outside_num_G = outside_num_G + data_R[c_indx] - data_G[c_indx + h_diff];
                        outside_num_B = outside_num_B + data_R[c_indx] - data_B[c_indx + h_diff];
                        inside_num_R = inside_num_R - data_R[c_indx] + data_R[c_indx + h_diff];
                        inside_num_G = inside_num_G - data_G[c_indx] + data_G[c_indx + h_diff];
                        inside_num_B = inside_num_B - data_B[c_indx] + data_B[c_indx + h_diff];
                    }
                    int lim = w - i * 8;
                    for (int n = 0; n < lim; n++) {
                        int c_indx = i + row_v_num * (up_left_c_y - 1);
                        int h_diff = row_v_num * h;
                        outside_num_R[n] = outside_num_R[n] + data_R[c_indx][n] - data_R[c_indx + h_diff][n];
                        outside_num_G[n] = outside_num_G[n] + data_R[c_indx][n] - data_G[c_indx + h_diff][n];
                        outside_num_B[n] = outside_num_B[n] + data_R[c_indx][n] - data_B[c_indx + h_diff][n];
                        inside_num_R[n] = inside_num_R[n] - data_R[c_indx][n] + data_R[c_indx + h_diff][n];
                        inside_num_G[n] = inside_num_G[n] - data_G[c_indx][n] + data_G[c_indx + h_diff][n];
                        inside_num_B[n] = inside_num_B[n] - data_B[c_indx][n] + data_B[c_indx + h_diff][n];
                    }
                }

                // Save values for movement in x direction
                double8_t tmp_outside_num_R = outside_num_R;
                double8_t tmp_outside_num_G = outside_num_G;
                double8_t tmp_outside_num_B = outside_num_B;
                double8_t tmp_inside_num_R = inside_num_R;
                double8_t tmp_inside_num_G = inside_num_G;
                double8_t tmp_inside_num_B = inside_num_B;

                for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {
                    
                    int low_right_c_x = w + up_left_c_x;

                    if (up_left_c_x != 0) {
                        // Go through values that are left outside the new region our that are added to the new region after movement in x direction:
                        int row_v = (up_left_c_x - 1)/8;
                        for (int j = up_left_c_y; j < low_right_c_y; j++) {
                            int c_indx = row_v + row_v_num * j;
                            int w_diff = w/8;
                            int n = (up_left_c_x - 1) - row_v*8;
                            std::cout << "n: " << n << std::endl;
                            std::cout << "row_v_num: " << row_v_num << std::endl;
                            std::cout << "c_index: " << c_indx << std::endl;
                            std::cout << "row_v_num * ny: " << row_v_num * ny << std::endl;
                            std::cout << "w_diff: " << w_diff << std::endl;
                            tmp_outside_num_R[n] = tmp_outside_num_R[n] + data_R[c_indx][n] - data_R[c_indx + w_diff][n];
                            tmp_outside_num_G[n] = tmp_outside_num_G[n] + data_G[c_indx][n] - data_G[c_indx + w_diff][n];
                            tmp_outside_num_B[n] = tmp_outside_num_B[n] + data_B[c_indx][n] - data_B[c_indx + w_diff][n];
                            tmp_inside_num_R[n] = tmp_inside_num_R[n] - data_R[c_indx][n] + data_R[c_indx + w_diff][n];
                            tmp_inside_num_G[n] = tmp_inside_num_G[n] - data_G[c_indx][n] + data_G[c_indx + w_diff][n];
                            tmp_inside_num_B[n] = tmp_inside_num_B[n] - data_B[c_indx][n] + data_B[c_indx + w_diff][n];
                        } 
                    }

                    // Calculating cost:
                    double8_t inside_sqr_num = {0,0,0,0,0,0,0,0};
                    double8_t outside_sqr_num = {0,0,0,0,0,0,0,0};
                    double8_t inside_denumerator = {double(h * w), double(h * w), double(h * w), double(h * w), double(h * w), double(h * w), double(h * w), double(h * w)};
                    double8_t outside_denumerator = double8_t {double(nx * ny), double(nx * ny), double(nx * ny), double(nx * ny), double(nx * ny), double(nx * ny), double(nx * ny), double(nx * ny)} - inside_denumerator;
                    inside_sqr_num = inside_sqr_num + tmp_inside_num_R * tmp_inside_num_R;
                    inside_sqr_num = inside_sqr_num + tmp_inside_num_G * tmp_inside_num_G;
                    inside_sqr_num = inside_sqr_num + tmp_inside_num_B * tmp_inside_num_B;
                    outside_sqr_num = outside_sqr_num + tmp_outside_num_R * tmp_outside_num_R;
                    outside_sqr_num = outside_sqr_num + tmp_outside_num_G * tmp_outside_num_G;
                    outside_sqr_num = outside_sqr_num + tmp_outside_num_B * tmp_outside_num_B;
                    double8_t cost = sqr_sum -(inside_sqr_num/inside_denumerator) -(outside_sqr_num/outside_denumerator);

                    double final_cost = cost[0] + cost[1] + cost[2] + cost[3] + cost[4] + cost[5] + cost[6] + cost[7];

                    struct Result res;

                    res.x0 = up_left_c_x;
                    res.x1 = low_right_c_x;
                    res.y0 = up_left_c_y;
                    res.y1 = low_right_c_y;
                    res.inner[0] = (tmp_inside_num_R[0] + tmp_inside_num_R[1] + tmp_inside_num_R[2] + tmp_inside_num_R[3] + tmp_inside_num_R[4] + tmp_inside_num_R[5] + tmp_inside_num_R[6] + tmp_inside_num_R[7])/double(h * w);
                    res.inner[1] = (tmp_inside_num_G[0] + tmp_inside_num_G[1] + tmp_inside_num_G[2] + tmp_inside_num_G[3] + tmp_inside_num_G[4] + tmp_inside_num_G[5] + tmp_inside_num_G[6] + tmp_inside_num_G[7])/double(h * w);
                    res.inner[2] = (tmp_inside_num_B[0] + tmp_inside_num_B[1] + tmp_inside_num_B[2] + tmp_inside_num_B[3] + tmp_inside_num_B[4] + tmp_inside_num_B[5] + tmp_inside_num_B[6] + tmp_inside_num_B[7])/double(h * w);
                    res.outer[0] = (tmp_outside_num_R[0] + tmp_outside_num_R[1] + tmp_outside_num_R[2] + tmp_outside_num_R[3] + tmp_outside_num_R[4] + tmp_outside_num_R[5] + tmp_outside_num_R[6] + tmp_outside_num_R[7])/double(nx * ny);
                    res.outer[1] = (tmp_outside_num_G[0] + tmp_outside_num_G[1] + tmp_outside_num_G[2] + tmp_outside_num_G[3] + tmp_outside_num_G[4] + tmp_outside_num_G[5] + tmp_outside_num_G[6] + tmp_outside_num_G[7])/double(nx * ny);
                    res.outer[2] = (tmp_outside_num_B[0] + tmp_outside_num_B[1] + tmp_outside_num_B[2] + tmp_outside_num_B[3] + tmp_outside_num_B[4] + tmp_outside_num_B[5] + tmp_outside_num_B[6] + tmp_outside_num_B[7])/double(nx * ny);

                    std::pair<double,Result> new_pair = std::pair<double,Result> (final_cost, res);
                    tmp_results_2[tmp_results_2_it] = new_pair;
                    tmp_results_2_it++;
                }
            }


            std::pair<double,Result> new_pair = dummy;
            double minimum_cost = DBL_MAX;

            for(int k = 0; k < tmp_results_2_it; k++) {
                std::pair<double,Result> tmp = tmp_results_2[k];
                if (tmp.first < minimum_cost) {
                    minimum_cost = tmp.first;
                    new_pair.first = tmp.first;
                    new_pair.second = tmp.second;
                }
            }

            tmp_results[h] = new_pair;

        }


        std::pair<double,Result> new_pair = dummy;
        double minimum_cost = DBL_MAX;

        for(auto it = tmp_results.begin(); it!=tmp_results.end(); it++) {
            std::pair<double,Result> tmp = *it;
            if (tmp.first < minimum_cost) {
                minimum_cost = tmp.first;
                new_pair.first = tmp.first;
                new_pair.second = tmp.second;
            }
        }

        results[w] = new_pair;

    }


    struct Result res;
    double minimum_cost = DBL_MAX;

    for(auto it = results.begin(); it!=results.end(); it++) {
        std::pair<double,Result> tmp = *it;
        if (tmp.first < minimum_cost) {
            minimum_cost = tmp.first;
            res = tmp.second;
        }
    }

    return res;
}
*/