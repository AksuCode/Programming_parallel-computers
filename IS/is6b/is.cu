#include <cfloat>
#include <vector>
#include <cmath>

struct Result {
    int y0;
    int x0;
    int y1;
    int x1;
    float outer[3];
    float inner[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
*/
Result segment(int ny, int nx, const float *data) {

    // Sum of all pixel color component values stored in vector (R,G,B).
    std::vector<double> color_sum_v = std::vector<double> (3,0);
    double csqr_sum = 0;    // Squared sum for all pixels and colors
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int c_indx = 3 * i + 3 * nx * j;
            for (int n = 0; n < 3; n++) {
                color_sum_v[n] = color_sum_v[n] + data[n + c_indx];
                csqr_sum = csqr_sum + double(data[n + c_indx]) * double(data[n + c_indx]);
            }
        }
    }

    struct Result res;
    double minimum_cost = DBL_MAX;
    
    // Create rectangle segments of varying size.
    for (int w = 1; w <= nx; w++) {
        for (int h = 1; h <= ny; h++) {

            // Inside region
            std::vector<double> inside_numerator = std::vector<double> (3,0);   // Inside segment (R, G, B)
            for (int j = 0; j < h; j++) {
                for (int i = 0; i < w; i++) {
                    int c_indx = 3 * i + 3 * nx * j;
                    for (int n = 0; n < 3; n++) {
                        inside_numerator[n] = inside_numerator[n] + data[n + c_indx];
                    }
                }
            }

            // Outside region
            std::vector<double> outside_numerator = std::vector<double> (3,0);
            for (int n = 0; n < 3; n++) {
                outside_numerator[n] = color_sum_v[n] - inside_numerator[n];
            }

            // The rectangle is moved. Inside and outside region numerators and sqr sums are recalculated
            for (int up_left_c_y = 0; h + up_left_c_y <= ny; up_left_c_y++) {

                int low_right_c_y = h + up_left_c_y;   // Lower right corner y pos

                if (up_left_c_y != 0) {
                    // Go through values that are left outside the new region our that are added to the new region after movement in y direction:
                    for (int i = 0; i < w; i++) {
                        int c_indx = 3 * i + 3 * nx * (up_left_c_y - 1);
                        int h_diff = 3 * nx * h;
                        for (int n = 0; n < 3; n++) {
                            outside_numerator[n] = outside_numerator[n] + data[n + c_indx] - data[n + c_indx + h_diff];
                            inside_numerator[n] = inside_numerator[n] - data[n + c_indx] + data[n + c_indx + h_diff];
                        }
                    }
                }

                // Save values for movement in x direction
                std::vector<double> tmp_outside_numerator = outside_numerator;
                std::vector<double> tmp_inside_numerator = inside_numerator;

                for (int up_left_c_x = 0; w + up_left_c_x <= nx; up_left_c_x++) {
                    
                    int low_right_c_x = w + up_left_c_x;

                    if (up_left_c_x != 0) {
                        // Go through values that are left outside the new region our that are added to the new region after movement in x direction:
                        for (int j = up_left_c_y; j < low_right_c_y; j++) {
                            int c_indx = 3 * (up_left_c_x - 1) + 3 * nx * j;
                            int w_diff = 3 * w;
                            for (int n = 0; n < 3; n++) {
                                tmp_outside_numerator[n] = tmp_outside_numerator[n] + data[n+ c_indx] - data[n + c_indx + w_diff];
                                tmp_inside_numerator[n] = tmp_inside_numerator[n] - data[n+ c_indx] + data[n + c_indx + w_diff];
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