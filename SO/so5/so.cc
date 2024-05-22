#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <climits>

#include <omp.h>

#include <iostream>

typedef unsigned long long data_t;

int qSort(int start, int end, data_t * data) {

    data_t pivot = data[end - 1];
    /*
    if (end - start > 10) {
        int pivot_indx = start + (end-start)/2;
        std::nth_element(data + pivot_indx - 4, data + pivot_indx, data + pivot_indx + 5);
        pivot = data[pivot_indx];
        data[pivot_indx] = data[end - 1];
        data[end - 1] = pivot;
    }
    */

    int l = start - 1;
    int r = (end - 1) + 1;
    while (true) {

        do {
            l++;
        } while (data[l] < pivot);

        do {
            r--;
        } while (data[r] >= pivot && r > l);

        if (l >= r) break;

        data_t tmp = data[l];
        data[l] = data[r];
        data[r] = tmp;

    }

    data_t tmp = data[end - 1];
    data[end - 1] = data[l];
    data[l] = tmp;

    return l;
}


void quickSort(int start, int end, data_t * data, int depth) {

    if (depth < 0) {
        std::sort(data + start, data + end);
        return;
    }

    if (start >= end) return;

    int new_div_point = qSort(start, end, data);
    
    #pragma omp task
    quickSort(start, new_div_point, data, depth - 1);

    #pragma omp task
    quickSort(new_div_point + 1, end, data, depth - 1);

}



void psort(int n, data_t *data) {

    if (n < 2) return;

    int recursion_depth = 8;
    #pragma omp parallel
    #pragma omp single
    quickSort(0, n, data, recursion_depth);

}

/*
int qSort(int start, int end, data_t * data) {

    data_t pivot = data[end - 1];
    if (end - start > 10) {
        int pivot_indx = start + (end-start)/2;
        std::nth_element(data + pivot_indx - 4, data + pivot_indx, data + pivot_indx + 5);
        pivot = data[pivot_indx];
        data[pivot_indx] = data[end - 1];
        data[end - 1] = pivot;
    }

    int l = start;
    for (int k = start; k < end - 1; k++) {
        if (data[k] <= pivot) {
            data_t tmp = data[k];
            data[k] = data[l];
            data[l] = tmp;
            l++;
        }
    }

    data_t tmp = data[l];
    data[l] = data[end - 1];
    data[end - 1] = tmp;

    return l;
}


void quickSort(int start, int end, data_t * data, int depth) {

    if (depth < 0) {
        std::sort(data + start, data + end);
        return;
    }

    if (start >= end) return;

    int new_div_point = qSort(start, end, data);
    
    #pragma omp task
    quickSort(start, new_div_point, data, depth - 1);

    #pragma omp task
    quickSort(new_div_point + 1, end, data, depth - 1);

}



void psort(int n, data_t *data) {

    if (n < 2) return;

    int recursion_depth = 8;
    #pragma omp parallel
    #pragma omp single
    quickSort(0, n, data, recursion_depth);

}
*/


/*
int qSort(int start, int end, data_t * data, data_t * aux) {

    data_t pivot_val = (data[start] + data[start + (end - start)/2] + data[end - 1])/3;

    int possible_div_point = start;
    data_t smallest_lefty = ULLONG_MAX;

    int right_it = end - 1;
    int left_it = start;
    for (int k = start; k < end; k++) {
        if (data[k] < pivot_val) {
            aux[left_it] = data[k];
            left_it++;
        } else {
            aux[right_it] = data[k];
            if (aux[right_it] < smallest_lefty) {
                smallest_lefty = aux[right_it];
                possible_div_point = right_it;
            }
            right_it--;
        }
    }

    right_it++;

    aux[possible_div_point] = aux[right_it];
    aux[right_it] = smallest_lefty;

    std::copy(aux + start, aux + end, data + start);

    return right_it;
}



void quickSort(int start, int end, data_t * data, data_t * aux, int depth) {

    if (depth < 0) {
        std::sort(data + start, data + end);
        return;
    }
    if (start >= end - 1) return;

    int new_div_point = qSort(start, end, data, aux);
    
    #pragma omp task
    quickSort(start, new_div_point, data, aux, depth - 1);

    #pragma omp task
    quickSort(new_div_point + 1, end, data, aux, depth - 1);

}



void psort(int n, data_t *data) {

    if (n < 2) return;

    int recursion_depth = 5;
    std::vector<data_t> aux(n);
    #pragma omp parallel
    #pragma omp single
    quickSort(0, n, data, &(aux[0]), recursion_depth);

}
*/