#include "matrix_pencil_method.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "fixed_point.h"
#include "matrix.h"
#include "svd_3_v2.h"

fixed_point_t sampling_frequency;

fixed_point_t L_parameter;  // pencil parameter

size_t M_parameter;  // number of eigenvalues

matrix_t data_matrix;

matrix_t data_matrix_v;

matrix_t computing_matrix;

// Print a matrix. Useful for debugging purposes.
void matrix_print(const matrix_t* matrix) {
    fixed_point_t entry = fixed_point_init(0);
    for (size_t i = 0; i < matrix->rows; ++i) {
        for (size_t j = 0; j < matrix->cols; ++j) {
            matrix_get(matrix, i, j, &entry);
            printf("%lld  ", entry);
        }
        printf("\n");
    }
    printf("\n");
}

// Determine the more positive eigenvalue of a 2x2 matrix.
// Return whether the computation was successful.
bool eigenvalues(const matrix_t* matrix, fixed_point_t* eigenvalue) {
    // Applying the quadratic formula to find the eigenvalues.
    fixed_point_t a = fixed_point_init(0);
    fixed_point_t b = fixed_point_init(0);
    fixed_point_t c = fixed_point_init(0);
    fixed_point_t d = fixed_point_init(0);

    matrix_get(matrix, 0, 0, &a);
    matrix_get(matrix, 0, 1, &b);
    matrix_get(matrix, 1, 0, &c);
    matrix_get(matrix, 1, 1, &d);

    fixed_point_t determinant = fixed_point_subtract(
        fixed_point_multiply(a, d), fixed_point_multiply(b, c));
    fixed_point_t n = fixed_point_add(a, d);
    fixed_point_t discriminant = fixed_point_subtract(
        fixed_point_square(n),
        fixed_point_multiply(fixed_point_init(4), determinant));

    if (discriminant < fixed_point_init(0)) {
        *eigenvalue = fixed_point_init(1);
        printf("Imaginary Eigenvalue\n");
    }

    fixed_point_t eigenvalue_1 = fixed_point_divide(
        fixed_point_add(n, fixed_point_square_root(discriminant)),
        fixed_point_init(2));
    fixed_point_t eigenvalue_2 = fixed_point_divide(
        fixed_point_add(n, fixed_point_square_root(discriminant)),
        fixed_point_init(2));

    // of the two eigenvalues, pick the smaller one that is not negative.
    if (eigenvalue_1 < fixed_point_init(0) || eigenvalue_2 < eigenvalue_1) {
        if (eigenvalue_2 < fixed_point_init(0)) {
            printf("Two Negative Eigenvalues\n");
        }
        *eigenvalue = eigenvalue_2;
    } else {
        *eigenvalue = eigenvalue_1;
    }

    return true;
}

// Form the (Y^T)Y matrix. We need to find the right singular vectors of Y.
bool form_data_matrix(int16_t data[], const size_t length, matrix_t* result) {
    fixed_point_t sum_of_squares = fixed_point_init(0);
    fixed_point_t sum_of_offset_1 = fixed_point_init(0);
    fixed_point_t sum_of_offset_2 = fixed_point_init(0);

    int16_t minimum_value_int = data[0];
    for (size_t i = 0; i < length; ++i) {
        if (data[i] < minimum_value_int) {
            minimum_value_int = data[i];
        }
    }

    fixed_point_t minimum_value = fixed_point_init(minimum_value_int);

    for (size_t i = 0; i < length; ++i) {
        sum_of_squares = fixed_point_add(
            sum_of_squares,
            fixed_point_multiply(fixed_point_subtract(fixed_point_init(data[i]),
                                                      minimum_value) >>
                                     6,
                                 fixed_point_subtract(fixed_point_init(data[i]),
                                                      minimum_value) >>
                                     6));
        if (i + 1 < length) {
            sum_of_offset_1 = fixed_point_add(
                sum_of_offset_1,
                fixed_point_multiply(
                    fixed_point_subtract(fixed_point_init(data[i]),
                                         minimum_value) >>
                        6,
                    fixed_point_subtract(fixed_point_init(data[i + 1]),
                                         minimum_value) >>
                        6));
        }
        if (i + 2 < length) {
            sum_of_offset_2 = fixed_point_add(
                sum_of_offset_2,
                fixed_point_multiply(
                    fixed_point_subtract(fixed_point_init(data[i]),
                                         minimum_value) >>
                        6,
                    fixed_point_subtract(fixed_point_init(data[i + 2]),
                                         minimum_value) >>
                        6));
        }
    }

    // form (Y^T)Y
    matrix_set(
        &data_matrix, 0, 0,
        fixed_point_subtract(
            fixed_point_subtract(
                sum_of_squares,
                fixed_point_multiply(
                    fixed_point_subtract(fixed_point_init(data[length - 2]),
                                         minimum_value) >>
                        6,
                    fixed_point_subtract(fixed_point_init(data[length - 2]),
                                         minimum_value) >>
                        6)),
            fixed_point_multiply(
                fixed_point_subtract(fixed_point_init(data[length - 1]),
                                     minimum_value) >>
                    6,
                fixed_point_subtract(fixed_point_init(data[length - 1]),
                                     minimum_value) >>
                    6)));
    matrix_set(&data_matrix, 0, 1,
               fixed_point_subtract(
                   sum_of_offset_1,
                   fixed_point_multiply(
                       fixed_point_subtract(fixed_point_init(data[length - 2]),
                                            minimum_value) >>
                           6,
                       fixed_point_subtract(fixed_point_init(data[length - 1]),
                                            minimum_value) >>
                           6)));
    matrix_set(&data_matrix, 0, 2, sum_of_offset_2);
    matrix_set(&data_matrix, 1, 0,
               fixed_point_subtract(
                   sum_of_offset_1,
                   fixed_point_multiply(
                       fixed_point_subtract(fixed_point_init(data[length - 2]),
                                            minimum_value) >>
                           6,
                       fixed_point_subtract(fixed_point_init(data[length - 1]),
                                            minimum_value) >>
                           6)));
    matrix_set(&data_matrix, 1, 1,
               fixed_point_subtract(
                   fixed_point_subtract(
                       sum_of_squares,
                       fixed_point_multiply(
                           fixed_point_subtract(fixed_point_init(data[0]),
                                                minimum_value) >>
                               6,
                           fixed_point_subtract(fixed_point_init(data[0]),
                                                minimum_value) >>
                               6)),
                   fixed_point_multiply(
                       fixed_point_subtract(fixed_point_init(data[length - 1]),
                                            minimum_value) >>
                           6,
                       fixed_point_subtract(fixed_point_init(data[length - 1]),
                                            minimum_value) >>
                           6)));
    matrix_set(
        &data_matrix, 1, 2,
        fixed_point_subtract(
            sum_of_offset_1,
            fixed_point_multiply(fixed_point_subtract(fixed_point_init(data[0]),
                                                      minimum_value) >>
                                     6,
                                 fixed_point_subtract(fixed_point_init(data[1]),
                                                      minimum_value) >>
                                     6)));
    matrix_set(&data_matrix, 2, 0, sum_of_offset_2);
    matrix_set(
        &data_matrix, 2, 1,
        fixed_point_subtract(
            sum_of_offset_1,
            fixed_point_multiply(fixed_point_subtract(fixed_point_init(data[0]),
                                                      minimum_value) >>
                                     6,
                                 fixed_point_subtract(fixed_point_init(data[1]),
                                                      minimum_value) >>
                                     6)));
    matrix_set(
        &data_matrix, 2, 2,
        fixed_point_subtract(
            fixed_point_subtract(
                sum_of_squares,
                fixed_point_multiply(
                    fixed_point_subtract(fixed_point_init(data[0]),
                                         minimum_value) >>
                        6,
                    fixed_point_subtract(fixed_point_init(data[0]),
                                         minimum_value) >>
                        6)),
            fixed_point_multiply(fixed_point_subtract(fixed_point_init(data[1]),
                                                      minimum_value) >>
                                     6,
                                 fixed_point_subtract(fixed_point_init(data[1]),
                                                      minimum_value) >>
                                     6)));

    return true;
}

// Form (((V1')^H)^t)((V2')^H) matrix. We need to find the eigenvalues of this
// matrix.
bool form_computing_matrix(void) {
    fixed_point_t a = fixed_point_init(0);
    fixed_point_t b = fixed_point_init(0);
    fixed_point_t c = fixed_point_init(0);
    fixed_point_t d = fixed_point_init(0);
    fixed_point_t e = fixed_point_init(0);
    fixed_point_t f = fixed_point_init(0);

    matrix_get(&data_matrix_v, 0, 0, &a);
    matrix_get(&data_matrix_v, 0, 1, &b);
    matrix_get(&data_matrix_v, 1, 0, &c);
    matrix_get(&data_matrix_v, 1, 1, &d);
    matrix_get(&data_matrix_v, 2, 0, &e);
    matrix_get(&data_matrix_v, 2, 1, &f);

    matrix_set(&computing_matrix, 0, 0, fixed_point_init(0));
    matrix_set(
        &computing_matrix, 0, 1,
        fixed_point_divide(fixed_point_subtract(fixed_point_multiply(d, e),
                                                fixed_point_multiply(c, f)),
                           fixed_point_subtract(fixed_point_multiply(a, d),
                                                fixed_point_multiply(b, c))));
    matrix_set(&computing_matrix, 1, 0, fixed_point_init(1));
    matrix_set(
        &computing_matrix, 1, 1,
        fixed_point_divide(fixed_point_subtract(fixed_point_multiply(a, f),
                                                fixed_point_multiply(b, e)),
                           fixed_point_subtract(fixed_point_multiply(a, d),
                                                fixed_point_multiply(b, c))));

    return true;
}

bool matrix_pencil_method_init(const fixed_point_t signal_sampling_frequency,
                               int16_t data[], const size_t length) {
    sampling_frequency = signal_sampling_frequency;
    L_parameter = fixed_point_init(2);
    M_parameter = 2;

    if (!matrix_init(&data_matrix, 3, 3)) {
        return false;
    }

    form_data_matrix(data, (3 * length) / 5, &data_matrix);

    if (!svd_init(&data_matrix)) {
        return false;
    }

    if (!matrix_init(&data_matrix_v, 3, 3)) {
        return false;
    }

    svd_calculate_v(&data_matrix_v);
    free(&data_matrix);

    if (!matrix_init(&computing_matrix, 2, 2)) {
        return false;
    }

    form_computing_matrix();

    return true;
}

// Compute e^y and return the result.
fixed_point_t exp(const fixed_point_t y) {
    fixed_point_t sum = fixed_point_init(1);
    fixed_point_t term = fixed_point_init(1);
    for (int i = 1; i < 20; ++i) {
        term = fixed_point_multiply(term,
                                    fixed_point_divide(y, fixed_point_init(i)));
        sum = fixed_point_add(sum, term);
    }
    return sum;
}

// Determine the natural log of a fixed_point number x.
// Return whether the computation was successful.
bool ln(const fixed_point_t x, fixed_point_t* result) {
    if (x <= fixed_point_init(0)) {
        return false;
    }

    fixed_point_t y = fixed_point_subtract(x, fixed_point_init(1));
    fixed_point_t term = fixed_point_init(1);

    while (fixed_point_absolute_value(fixed_point_multiply(
               term, fixed_point_init(2048))) > fixed_point_init(1)) {
        term = fixed_point_subtract(fixed_point_divide(x, exp(y)),
                                    fixed_point_init(1));
        y = fixed_point_add(y, term);
    }

    *result = y;
    return true;
}

bool matrix_pencil_method_get_time_constant(fixed_point_t* time_constant) {
    fixed_point_t eigenvalue = fixed_point_init(0);
    fixed_point_t n = fixed_point_init(0);

    if (!eigenvalues(&computing_matrix, &eigenvalue)) {
        return false;
    }

    if (!ln(eigenvalue, &n)) {
        return false;
    }

    *time_constant = fixed_point_absolute_value(fixed_point_divide(
        fixed_point_init(1), fixed_point_multiply(sampling_frequency, n)));

    return true;
}
