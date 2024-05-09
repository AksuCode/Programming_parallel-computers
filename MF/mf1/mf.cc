#include <vector>
#include <algorithm>
#include <iostream>

/*
This is the function you need to implement. Quicn reference:
- input rows: 0 <= y < ny
- input columns: 0 <= x < nx
- element at row y and column x is stored in in[x + y*nx]
- for each pixel (x, y), store the median of the pixels (a, b) which satisfy
  max(x-hx, 0) <= a < min(x+hx+1, nx), max(y-hy, 0) <= b < min(y+hy+1, ny)
  in out[x + y*nx].
*/
void mf(int ny, int nx, int hy, int hx, const float *in, float *out) {

  for (int y = 0; y < ny; y++) {
    for (int x = 0; x < nx; x++) {

      std::vector<float> aux((2*hx+1) * (2*hy+1), 0);

      int j = (y - hy);
      if (j < 0) j = 0;
      int n = 0;
      for (; j < ny && j <= (y + hy); j++) {
        int i = (x - hx);
        if (i < 0) i = 0;
        for(; i < nx && i <= (x + hx); i++){
          aux[n] = in[i + nx*j];
          n++;
        }
      }


      int halved = n/2;
      std::nth_element(aux.begin(), aux.begin() + halved, aux.begin() + n);
      if (n % 2 == 0) {
        float step = aux[halved];
        std::nth_element(aux.begin(), aux.begin() + halved - 1, aux.begin() + n);
        out[x + nx*y] = (aux[halved - 1] + step)/2;  // - 1 because vector starts with index 0
      } else {
        out[x + nx*y] = aux[halved];
      }

    }
  }

}
