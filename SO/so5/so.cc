#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <climits>

#include <omp.h>

#include <iostream>

typedef unsigned long long data_t;

int qSort(int start, int end, data_t * data) {

    data_t pivot = data[start + (end-start)/2];

    /*
    // Median value of first, mid and last
    int pivot_indx = start;
    data_t pivot = data[start];

    if (data[start + (end-start)/2] > pivot) {
        if (data[end - 1] > data[start + (end-start)/2]) {
            pivot = data[start + (end-start)/2];
            //ivot_indx = start + (end-start)/2;
        } else {
            if (pivot < data[end - 1]) {
                pivot = data[end - 1];
                //pivot_indx = end - 1;
            }
        }
    } else {
        if (pivot > data[end - 1]) {
            pivot = data[end - 1];
            //pivot_indx = end - 1;
        }
    }
    */

    int l = start;
    for (; l < end; l++) {
        if (data[l] > pivot) break;
    }

    for (int r = l + 1; r < end; r++) {
        if (data[r] <= pivot) {
            data_t tmp = data[r];
            data[r] = data[l];
            data[l] = tmp;
            l++;
        }
    }

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
    quickSort(new_div_point, end, data, depth - 1);

}



void psort(int n, data_t *data) {

    if (n < 2) return;

    int recursion_depth = 5;
    #pragma omp parallel
    #pragma omp single
    quickSort(0, n, data, recursion_depth);

}


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