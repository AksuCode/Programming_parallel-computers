#include <algorithm>
#include <vector>

#include <omp.h>

#include <iostream>

typedef unsigned long long data_t;

/*
void mSort(int start, int mid, int end, data_t * data) {
    int arr_size = end - start;
    std::vector<data_t> aux(arr_size);
    int right = mid;
    int left = start;

    for (int k = 0; k < arr_size; k++) {

        if (right == end) {
            aux[k] = data[left];
            left++;
            continue;
        }
        
        if (left == mid) {
            aux[k] = data[right];
            right++;
            continue;
        }


        if (data[right] < data[left]) {
            aux[k] = data[right];
            right++;
        } else {
            aux[k] = data[left];
            left++;
        }
    }

    data_t * ptr = data + start;
    for (auto it = aux.begin(); it!=aux.end(); it++) {
        *ptr = *it;
        ptr++;
    }
}


void mergeSort(int start, int end, data_t * data, int depth) {

    if (depth <= 0) {
        std::sort(data + start, data + end);
        return;
    }

    int mid = start + (end - start)/2;
    if (mid <= start) return;
    
    mergeSort(start, mid, data, depth - 1);
    mergeSort(mid, end, data, depth - 1);

    mSort(start, mid, end, data);

}

void psort(int n, data_t *data) {

    int recursion_depth = 4;
    mergeSort(0, n, data, recursion_depth);

}
*/












/*
void mSort(data_t *start, data_t *mid, data_t *end) {
    data_t * right = mid;
    for (data_t * ptr = start; ptr != end && ptr != right && right != end; ptr++) {
        if (*right < *ptr) {
            data_t tmp =*ptr;
            *ptr = *right;
            *right = tmp;
            right++;
        }
    }
}
*/


/*
void psort(int n, data_t *data) {

    int thread_c = 20;
    int mid = n/thread_c;
    #pragma omp parallel for
    for (int i = 0; i < thread_c; i++) {
        std::sort(data + mid*i, data + mid*i + mid);
    }
    std::sort(data, data + n);
}
*/

/*
void psort(int n, data_t *data) {

    std::sort(data, data + n);
}
*/

/*
void mergeSort(int start, int end, data_t * data, int depth) {
    if (depth <= 0) return;

    int mid = start + (end - start)/2;
    if (mid <= start) return;

    #pragma omp task
    mergeSort(start, mid, data, depth - 1);
    #pragma omp task
    mergeSort(mid, end, data, depth - 1);

    #pragma omp critical
    {
        std::sort(data + start, data + end);
    }

}

void psort(int n, data_t *data) {

    int recursion_depth = 10;
    #pragma omp parallel
    #pragma omp single
    {
        mergeSort(0, n, data, recursion_depth);
    }
    std::sort(data, data + n);

}
*/


/*
void mergeSort(int start, int end, data_t * data) {
    int interval = end - start;
    int mid = start + interval/2;
    
    if (mid == start) {
        return;
    }

    if (start < end - 1) {
        mergeSort(start, mid, data);
        mergeSort(mid, end, data);
    }

    std::cout << "Before merge: " << std::endl;
    for(int k = start; k < end; k++) {
        std::cout << data[k] << std::endl;
    }
    std::cout << std::endl;

    if (data[mid] < data[start]) {
        int halved_int = interval/2;
        if (interval % 2 == 0)  {
            for (int k = 0; k < halved_int; k++) {
                data_t tmp = data[start + k];
                data[start + k] = data[mid + k];
                data[mid + k] = tmp;
            }
        } else {
            data_t tmp_1 = data[mid];
            for (int k = 0; k < halved_int; k++) {
                data_t tmp_2 = data[mid + k + 1];
                data[mid + k + 1] = data[start + k];
                data[start + k] = tmp_1;
                tmp_1 = tmp_2;
            }
            data[mid] = tmp_1;
        }
    }

    std::cout << "After merge: " << std::endl;
    for(int k = start; k < end; k++) {
        std::cout << data[k] << std::endl;
    }
    std::cout << std::endl;

    std::sort(data + start, data + end);


    std::cout << "After sort: " << std::endl;
    for(int k = start; k < end; k++) {
        std::cout << data[k] << std::endl;
    }
    std::cout << std::endl;


}

void psort(int n, data_t *data) {

    if (n == 1) return;

    mergeSort(0, n, data);

    std::sort(data, data + n);

}
*/

void mSort(int start, int mid, int end, data_t * data) {
    int arr_size = end - start;
    std::vector<data_t> aux(arr_size);
    int right = mid;
    int left = start;

    for (int k = start; k < mid; k++) {
        std::cout << data[k] << std::endl;
    }
    std::cout << std::endl;

    for (int k = mid; k < end; k++) {
        std::cout << data[k] << std::endl;
    }
    std::cout << std::endl;

    for (int k = 0; k < arr_size; k++) {

        if (right == end) {
            aux[k] = data[left];
            left++;
            continue;
        }
        
        if (left == mid) {
            aux[k] = data[right];
            right++;
            continue;
        }


        if (data[right] < data[left]) {
            aux[k] = data[right];
            right++;
        } else {
            aux[k] = data[left];
            left++;
        }
    }

    data_t * ptr = data + start;
    for (auto it = aux.begin(); it!=aux.end(); it++) {
        *ptr = *it;
        ptr++;
    }
}

void psort(int n, data_t *data) {

    const int depth = 1;
    for (int k = depth; k <= n; k=k*2) {
        int i = 0;
        for (; i < (n/k) - 1; i++) {
            int start = k*i;
            int end = k*(i+1);

            if (k == depth) {
                std::sort(data + start, data +  start + (end - start)/2);
                std::sort(data + start + (end - start)/2, data + end);
            }
            else mSort(start, start + (end - start)/2, end, data);
        }
        int start = k*i;
        int end = n;

        if (k == depth) {
            std::sort(data + start, data + start + (end - start)/2);
            std::sort(data + start + (end - start)/2, data + end);
        }
        else mSort(start, start + (end - start)/2, end, data);
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