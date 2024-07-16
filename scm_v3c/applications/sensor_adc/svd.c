
// ONLY FOR 2X2 MATRICES

#include "svd.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "matrix.h"

matrix_t svd_computing_matrix;

bool svd_init(const matrix_t* matrix) {
    if (matrix->cols != 2 || matrix->rows != 2) {
        return false;
    }
    matrix_copy(matrix, &svd_computing_matrix);
    return true;
}

// Calculates the square root of s with relative tolerance 0.0001.
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Heron's_method
matrix_type_t herons_method(matrix_type_t s) {
    matrix_type_t guess = s;
    while ((guess * guess) - s > 0.0001 * s) {
        guess = (guess + (s / guess)) / 2;
    }
    return guess;
}

// Calculates the square of a number n.
matrix_type_t square(matrix_type_t n) {
    if (n < 0) {
        matrix_type_t positive_n = -n;
        return positive_n * positive_n;
    }
    return n * n;
}

// Calculate the eigenvalues of a 2x2 matrix with entries a, b, c, and d.
void eigenvalues(matrix_type_t a, matrix_type_t b, matrix_type_t c,
                 matrix_type_t d, matrix_type_t* eigenvalue_1,
                 matrix_type_t* eigenvalue_2) {
    // Applying the quadratic formula to find the eigenvalues.
    *eigenvalue_1 =
        (a + d + herons_method(square(a + d) - (4 * (square(a) - square(b))))) /
        2;
    *eigenvalue_2 =
        (a + d - herons_method(square(a + d) - (4 * (square(a) - square(b))))) /
        2;
}

// Calculate the nullspace of a 2x2 matrix.
// Returns whether the calculation was successful.
// https://www.quora.com/How-do-you-find-the-null-space-of-a-2x2-matrix
void nullspace(matrix_type_t a, matrix_type_t b, matrix_type_t c,
               matrix_type_t d, matrix_type_t* nullspace_entry_1,
               matrix_type_t* nullspace_entry_2) {
    if ((square(a) - square(b)) != 0) {
        // If the determinant of the matrix does not equal 0, then it has a
        // trivial nullspace.
        *nullspace_entry_1 = 0;
        *nullspace_entry_2 = 0;
        return;
    }
    if (a != 0) {
        *nullspace_entry_1 = -b;
        *nullspace_entry_2 = a;
        return;
    }
    if (c != 0) {  // a = 0
        *nullspace_entry_1 = -d;
        *nullspace_entry_2 = c;
        return;
    }
    if (b == 0 && d == 0) {  // a = 0 and c = 0
        // This is the zero matrix case, and all vectors in R^2 are in the
        // nullspace, so just return [1, 1]^T (i.e. one vector).
        *nullspace_entry_1 = 1;
        *nullspace_entry_2 = 1;
        return;
    }
    *nullspace_entry_1 = 1;
    *nullspace_entry_2 = 0;
    return;
}

// Get the entries of a 2x2 matrix in the form [a b
//                                              c d]
void get_matrix_entries(const matrix_t* matrix, matrix_type_t* a,
                        matrix_type_t* b, matrix_type_t* c, matrix_type_t* d) {
    matrix_get(matrix, 0, 0, a);
    matrix_get(matrix, 0, 1, b);
    matrix_get(matrix, 1, 0, c);
    matrix_get(matrix, 1, 1, d);
}

void svd_calculate_v(matrix_t* result) {
    matrix_type_t a = 0;
    matrix_type_t b = 0;
    matrix_type_t c = 0;
    matrix_type_t d = 0;
    get_matrix_entries(&svd_computing_matrix, &a, &b, &c, &d);

    // Computing with the matrix (A^T)A.
    matrix_set(&svd_computing_matrix, 0, 0, (a * a) + (c * c));
    matrix_set(&svd_computing_matrix, 0, 1, (a * b) + (c * d));
    matrix_set(&svd_computing_matrix, 1, 0, (a * b) + (d * c));
    matrix_set(&svd_computing_matrix, 1, 1, (b * b) + (d * d));

    // Get the entries of (A^T)A
    get_matrix_entries(&svd_computing_matrix, &a, &b, &c, &d);

    // Calculate the eigenvalues of (A^T)A.
    matrix_type_t eigenvalue_1 = -1;
    matrix_type_t eigenvalue_2 = -1;
    if (square(a + d) - (4 * (square(a) - square(b))) >= 0) {
        eigenvalues(a, b, c, d, &eigenvalue_1, &eigenvalue_2);
    } else {
        eigenvalue_1 = 0;
        eigenvalue_2 = 0;
    }

    // Calculate the eigenvectors of (A^T)A and add them to the V matrix result.
    matrix_type_t eigenvector_entry1, eigenvector_entry2;
    nullspace(a - eigenvalue_1, b, c, d - eigenvalue_1, &eigenvector_entry1,
              &eigenvector_entry2);
    matrix_set(result, 0, 0, eigenvector_entry1);
    matrix_set(result, 1, 0, eigenvector_entry2);

    nullspace(a - eigenvalue_2, b, c, d - eigenvalue_2, &eigenvector_entry1,
              &eigenvector_entry2);
    matrix_set(result, 0, 1, eigenvector_entry1);
    matrix_set(result, 1, 1, eigenvector_entry2);
}
