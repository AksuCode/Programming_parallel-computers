#include <algorithm>
#include <vector>

typedef unsigned long long data_t;

void psort(int n, data_t *data) {

    for (int k = 1; k <= n; k=k*2) {
        #pragma omp parallel for
        for (int i = 0; i < (n/k) - 1; i++) {
            int start = k*i;
            int end = start + k;

            std::sort(data + start, data + end);
        }
        int start = k*((n/k) - 1);
        int end = n;

        std::sort(data + start, data + end);
    }
}

/* WORKING simple implementation.
Should work just like merge sort but no additional memory (aux arrays) used:
void psort(int n, data_t *data) {

    for (int k = 1; k <= n; k=k*2) {
        int i = 0;
        for (; i < (n/k) - 1; i++) {
            int start = k*i;
            int end = k*(i+1);

            std::sort(data + start, data + end);
        }
        int start = k*i;
        int end = n;

        std::sort(data + start, data + end);
    }
}
*/