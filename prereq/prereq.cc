struct Result {
    float avg[3];
};

/*
This is the function you need to implement. Quick reference:
- x coordinates: 0 <= x < nx
- y coordinates: 0 <= y < ny
- horizontal position: 0 <= x0 < x1 <= nx
- vertical position: 0 <= y0 < y1 <= ny
- color components: 0 <= c < 3
- input: data[c + 3 * x + 3 * nx * y]
- output: avg[c]
*/

/*
Result calculate(int ny, int nx, const float *data, int y0, int x0, int y1, int x1) {
    double red = 0;
    double green = 0;
    double blue = 0;

    for (int y = y0; y < y1; y++) {
        for (int x = x0; x < x1; x++) {
            int block_of_c = 3 * (x + nx * y);
            red = red + data[block_of_c];
            green = green + data[1 + block_of_c];
            blue = blue + data[2 + block_of_c];
        }
    }

    int area = (y1 - y0) * (x1 - x0);

    Result result{{float(red/area), float(green/area), float(blue/area)}};
    return result;
}
*/

Result calculate(int ny, int nx, const float *data, int y0, int x0, int y1, int x1) {
    double red = 0;
    double green = 0;
    double blue = 0;

    for (int i = 3 * (x0 + nx * y0); i <= 3 * ((x1-1) + nx * (y1-1)); i = i + 3) {

        if (i == x1) {
            i = i + 3 * nx;
            continue;
        }
        
        red = red + data[i];
        green = green + data[1 + i];
        blue = blue + data[2 + i];

    }

    int area = (y1 - y0) * (x1 - x0);

    Result result{{float(red/area), float(green/area), float(blue/area)}};
    return result;
}