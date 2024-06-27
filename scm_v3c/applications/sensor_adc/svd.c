
// ONLY FOR 2X2 MATRICES

#include "svd.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "matrix.h"

matrix_t mtx;

bool svd_init(matrix_t* m) {
    if (m->cols != 2 || m->rows != 2) {
        return false;
    }
    matrix_init(&mtx, 2, 2);
    matrix_copy(m, &mtx);
    return true;
}

// calculates the square root of s with tolerance 0.0001
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Heron's_method 
matrix_type_t herons_method(matrix_type_t s) {
    matrix_type_t guess = s;
    while ((guess * guess) - s > 0.0001 * s) {
        guess = (guess + (s / guess)) / 2;
    }
    return guess;
}

// calculates the square of a number n
matrix_type_t square(matrix_type_t n) {
    if (n < 0) {
        matrix_type_t positive_n = -n;
        return positive_n * positive_n;
    } else {
        return n * n;
    }
}

// calculate the eigenvalues of a 2x2 matrix with entries a, b, c, and d
void eigenvalues(matrix_type_t a, matrix_type_t b, 
                matrix_type_t c, matrix_type_t d,
                matrix_type_t* eigenvalue_1, matrix_type_t* eigenvalue_2) {
    *eigenvalue_1 = (a + d + herons_method(square(a + d) - (4 * (square(a) - square(b))))) / 2;
    *eigenvalue_2 = (a + d - herons_method(square(a + d) - (4 * (square(a) - square(b))))) / 2;
}

// calculate the nullspace of a 2x2 matrix
// https://www.quora.com/How-do-you-find-the-null-space-of-a-2x2-matrix 
void nullspace(matrix_type_t a, matrix_type_t b, 
                matrix_type_t c, matrix_type_t d, 
                matrix_type_t* nullspace_entry_1, 
                matrix_type_t* nullspace_entry_2) {
    if ((square(a) - square(b)) != 0) {
        *nullspace_entry_1 = 0; 
        *nullspace_entry_2 = 0;
    } else if (a == 0) {
        if (c == 0){
            if (b == 0 && d == 0) { // all vectors in R^2 are eigenvectors, so just return [1, 1]^T
                *nullspace_entry_1 = 1; 
                *nullspace_entry_2 = 1;
            } else {
                *nullspace_entry_1 = 1; 
                *nullspace_entry_2 = 0;
            }
        } else {
            *nullspace_entry_1 = -d; 
            *nullspace_entry_2 = c;
        }
    } else {
        *nullspace_entry_1 = -b; 
        *nullspace_entry_2 = a;
    }
}

void get_matrix_entries(const matrix_t* matrix, matrix_type_t* a, 
                        matrix_type_t* b, matrix_type_t* c, matrix_type_t* d) {
    matrix_get(matrix, 0, 0, a);
    matrix_get(matrix, 0, 1, b);
    matrix_get(matrix, 1, 0, c);
    matrix_get(matrix, 1, 1, d);
}

// calculates the SVD of a 2x2 matrix
void svd_calculate_v(matrix_t* result) {
    matrix_type_t a, b, c, d;
    get_matrix_entries(&mtx, &a, &b, &c, &d);

    // computing with the matrix (A^T)A
    matrix_set(&mtx, 0, 0, (a * a) + (c * c));
    matrix_set(&mtx, 0, 1, (a * b) + (c * d));
    matrix_set(&mtx, 1, 0, (a * b) + (d * c));
    matrix_set(&mtx, 1, 1, (b * b) + (d *d));

    // 2x2 matrix of the form [a b
    //                         c d]
    get_matrix_entries(&mtx, &a, &b, &c, &d);

    // calculate the eigenvalues
    matrix_type_t eigenvalue_1, eigenvalue_2;
    if (square(a + d) - (4 * (square(a) - square(b))) >= 0) {
        eigenvalues(a, b, c, d, &eigenvalue_1, &eigenvalue_2);
    } else {
        eigenvalue_1 = 0;
        eigenvalue_2 = 0;
    }

    // calculate the eigenvectors and add them to the v matrix result
    matrix_type_t eigenvector_entry1, eigenvector_entry2;
    nullspace(a - eigenvalue_1, b, c, d - eigenvalue_1, &eigenvector_entry1, &eigenvector_entry2);
    matrix_set(result, 0, 0, eigenvector_entry1);
    matrix_set(result, 1, 0, eigenvector_entry2);
    
    nullspace(a - eigenvalue_2, b, c, d - eigenvalue_2, &eigenvector_entry1, &eigenvector_entry2);    
    matrix_set(result, 0, 1, eigenvector_entry1);
    matrix_set(result, 1, 1, eigenvector_entry2);
}
