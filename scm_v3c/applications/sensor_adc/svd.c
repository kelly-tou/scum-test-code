
// ONLY FOR 2X2 MATRICES

#include "svd.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "fixed_point.h"
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
fixed_point_short_t herons_method(fixed_point_short_t s) {
    fixed_point_short_t guess = s;
    fixed_point_short_t tolerance = fixed_point_init(0.0001);
    while (fixed_point_subtract(fixed_point_multiply(guess, guess), s) >
           fixed_point_multiply(tolerance, s)) {
        guess = fixed_point_divide(
            fixed_point_add(guess, fixed_point_divide(s, guess)),
            fixed_point_init(2));
    }
    return guess;
}

// Calculates the square of a number n.
fixed_point_short_t square(fixed_point_short_t n) {
    if (n < fixed_point_init(0)) {
        fixed_point_short_t positive_n = -n;
        return fixed_point_multiply(positive_n, positive_n);
    }
    return fixed_point_multiply(n, n);
}

// Calculate the eigenvalues of a 2x2 matrix with entries a, b, c, and d.
void eigenvalues(fixed_point_short_t a, fixed_point_short_t b,
                 fixed_point_short_t c, fixed_point_short_t d,
                 fixed_point_short_t* eigenvalue_1,
                 fixed_point_short_t* eigenvalue_2) {
    // Applying the quadratic formula to find the eigenvalues.
    fixed_point_short_t n = fixed_point_add(a, d);
    *eigenvalue_1 = fixed_point_divide(
        fixed_point_add(
            n,
            herons_method(fixed_point_subtract(
                square(n), fixed_point_multiply(
                               fixed_point_init(4),
                               fixed_point_subtract(square(a), square(b)))))),
        fixed_point_init(2));
    *eigenvalue_2 = fixed_point_divide(
        fixed_point_subtract(
            n,
            herons_method(fixed_point_subtract(
                square(n), fixed_point_multiply(
                               fixed_point_init(4),
                               fixed_point_subtract(square(a), square(b)))))),
        fixed_point_init(2));
}

// Calculate the nullspace of a 2x2 matrix.
// Returns whether the calculation was successful.
// https://www.quora.com/How-do-you-find-the-null-space-of-a-2x2-matrix
void nullspace(fixed_point_short_t a, fixed_point_short_t b,
               fixed_point_short_t c, fixed_point_short_t d,
               fixed_point_short_t* nullspace_entry_1,
               fixed_point_short_t* nullspace_entry_2) {
    if (fixed_point_subtract(square(a), square(b)) != fixed_point_init(0)) {
        // If the determinant of the matrix does not equal 0, then it has a
        // trivial nullspace.
        *nullspace_entry_1 = fixed_point_init(0);
        *nullspace_entry_2 = fixed_point_init(0);
        return;
    }
    if (a != fixed_point_init(0)) {
        *nullspace_entry_1 = -b;
        *nullspace_entry_2 = a;
        return;
    }
    if (c != fixed_point_init(0)) {  // a = 0
        *nullspace_entry_1 = -d;
        *nullspace_entry_2 = c;
        return;
    }
    if (b == fixed_point_init(0) &&
        d == fixed_point_init(0)) {  // a = 0 and c = 0
        // This is the zero matrix case, and all vectors in R^2 are in the
        // nullspace, so just return [1, 1]^T (i.e. one vector).
        *nullspace_entry_1 = fixed_point_init(1);
        *nullspace_entry_2 = fixed_point_init(1);
        return;
    }
    *nullspace_entry_1 = fixed_point_init(1);
    *nullspace_entry_2 = fixed_point_init(0);
    return;
}

// Get the entries of a 2x2 matrix in the form [a b
//                                              c d]
void get_matrix_entries(const matrix_t* matrix, fixed_point_short_t* a,
                        fixed_point_short_t* b, fixed_point_short_t* c,
                        fixed_point_short_t* d) {
    matrix_get(matrix, 0, 0, a);
    matrix_get(matrix, 0, 1, b);
    matrix_get(matrix, 1, 0, c);
    matrix_get(matrix, 1, 1, d);
}

void svd_calculate_v(matrix_t* result) {
    fixed_point_short_t a = fixed_point_init(0);
    fixed_point_short_t b = fixed_point_init(0);
    fixed_point_short_t c = fixed_point_init(0);
    fixed_point_short_t d = fixed_point_init(0);
    get_matrix_entries(&svd_computing_matrix, &a, &b, &c, &d);

    // Computing with the matrix (A^T)A.
    matrix_set(&svd_computing_matrix, 0, 0,
               fixed_point_add(square(a), square(c)));
    matrix_set(&svd_computing_matrix, 0, 1,
               fixed_point_add(fixed_point_multiply(a, b),
                               fixed_point_multiply(c, d)));
    matrix_set(&svd_computing_matrix, 1, 0,
               fixed_point_add(fixed_point_multiply(a, b),
                               fixed_point_multiply(d, c)));
    matrix_set(&svd_computing_matrix, 1, 1,
               fixed_point_add(square(b), square(d)));

    // Get the entries of (A^T)A
    get_matrix_entries(&svd_computing_matrix, &a, &b, &c, &d);

    // Calculate the eigenvalues of (A^T)A.
    fixed_point_short_t eigenvalue_1 = -fixed_point_init(1);
    fixed_point_short_t eigenvalue_2 = -fixed_point_init(1);
    if (fixed_point_subtract(
            square(fixed_point_add(a, d)),
            fixed_point_multiply(fixed_point_init(4),
                                 fixed_point_subtract(square(a), square(b)))) >=
        fixed_point_init(0)) {
        eigenvalues(a, b, c, d, &eigenvalue_1, &eigenvalue_2);
    } else {
        eigenvalue_1 = fixed_point_init(0);
        eigenvalue_2 = fixed_point_init(0);
    }

    // Calculate the eigenvectors of (A^T)A and add them to the V matrix result.
    fixed_point_short_t eigenvector_entry1, eigenvector_entry2;
    nullspace(fixed_point_subtract(a, eigenvalue_1), b, c,
              fixed_point_subtract(d, eigenvalue_1), &eigenvector_entry1,
              &eigenvector_entry2);
    matrix_set(result, 0, 0, eigenvector_entry1);
    matrix_set(result, 1, 0, eigenvector_entry2);

    nullspace(fixed_point_subtract(a, eigenvalue_2), b, c,
              fixed_point_subtract(d, eigenvalue_2), &eigenvector_entry1,
              &eigenvector_entry2);
    matrix_set(result, 0, 1, eigenvector_entry1);
    matrix_set(result, 1, 1, eigenvector_entry2);
}
