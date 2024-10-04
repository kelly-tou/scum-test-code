#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "fixed_point.h"
#include "matrix.h"

matrix_t AtA;

bool svd_init(const matrix_t* matrix) {
    AtA = *matrix;
    return true;
}

bool cross_product(const matrix_t* matrix, matrix_t* result) {
    if (matrix->rows != 3) {
        return false;
    }
    fixed_point_t vector_1_entry_1 = fixed_point_init(0);
    fixed_point_t vector_1_entry_2 = fixed_point_init(0);
    fixed_point_t vector_1_entry_3 = fixed_point_init(0);
    fixed_point_t vector_2_entry_1 = fixed_point_init(0);
    fixed_point_t vector_2_entry_2 = fixed_point_init(0);
    fixed_point_t vector_2_entry_3 = fixed_point_init(0);

    matrix_get(matrix, 0, 0, &vector_1_entry_1);
    matrix_get(matrix, 1, 0, &vector_1_entry_2);
    matrix_get(matrix, 2, 0, &vector_1_entry_3);
    matrix_get(matrix, 0, 1, &vector_2_entry_1);
    matrix_get(matrix, 1, 1, &vector_2_entry_2);
    matrix_get(matrix, 2, 1, &vector_2_entry_3);

    matrix_set(result, 0, 0,
               fixed_point_subtract(
                   fixed_point_multiply(vector_1_entry_2, vector_2_entry_3),
                   fixed_point_multiply(vector_1_entry_3, vector_2_entry_2)));
    matrix_set(result, 1, 0,
               fixed_point_subtract(
                   fixed_point_multiply(vector_1_entry_3, vector_2_entry_1),
                   fixed_point_multiply(vector_1_entry_1, vector_2_entry_3)));
    matrix_set(result, 2, 0,
               fixed_point_subtract(
                   fixed_point_multiply(vector_1_entry_1, vector_2_entry_2),
                   fixed_point_multiply(vector_1_entry_2, vector_2_entry_1)));

    return true;
}

bool identity_subtract_outer_product(const matrix_t* vector, matrix_t* result) {
    if (vector->rows != 3 || vector->cols != 1 || result->rows != 3 ||
        result->cols != 3) {
        return false;
    }

    fixed_point_t vector_entry_1 = fixed_point_init(0);
    fixed_point_t vector_entry_2 = fixed_point_init(0);
    fixed_point_t vector_entry_3 = fixed_point_init(0);

    matrix_get(vector, 0, 0, &vector_entry_1);
    matrix_get(vector, 1, 0, &vector_entry_2);
    matrix_get(vector, 2, 0, &vector_entry_3);

    matrix_set(result, 0, 0,
               fixed_point_subtract(
                   fixed_point_init(1),
                   fixed_point_multiply(vector_entry_1, vector_entry_1)));
    matrix_set(result, 0, 1,
               -fixed_point_multiply(vector_entry_1, vector_entry_2));
    matrix_set(result, 0, 2,
               -fixed_point_multiply(vector_entry_1, vector_entry_3));
    matrix_set(result, 1, 0,
               -fixed_point_multiply(vector_entry_2, vector_entry_1));
    matrix_set(result, 1, 1,
               fixed_point_subtract(
                   fixed_point_init(1),
                   fixed_point_multiply(vector_entry_2, vector_entry_2)));
    matrix_set(result, 1, 2,
               -fixed_point_multiply(vector_entry_2, vector_entry_3));
    matrix_set(result, 2, 0,
               -fixed_point_multiply(vector_entry_3, vector_entry_1));
    matrix_set(result, 2, 1,
               -fixed_point_multiply(vector_entry_3, vector_entry_2));
    matrix_set(result, 2, 2,
               fixed_point_subtract(
                   fixed_point_init(1),
                   fixed_point_multiply(vector_entry_3, vector_entry_3)));

    return true;
}

void matrix_get_column(const matrix_t* matrix, const int8_t col,
                       matrix_t* result) {
    matrix_init(result, matrix->rows, 1);
    fixed_point_t entry = fixed_point_init(0);
    for (int8_t i = 0; i < matrix->rows; ++i) {
        matrix_get(matrix, i, col, &entry);
        matrix_set(result, i, 0, entry);
    }
}

fixed_point_t euclidean_norm(const matrix_t* vector) {
    fixed_point_t norm_squared = fixed_point_init(0);
    fixed_point_t entry = fixed_point_init(0);
    for (int8_t i = 0; i < vector->rows; ++i) {
        matrix_get(vector, i, 0, &entry);
        norm_squared = fixed_point_add(norm_squared, fixed_point_square(entry));
    }
    return fixed_point_square_root(norm_squared);
}

void normalize(matrix_t* vector) {
    fixed_point_t entry = fixed_point_init(0);
    matrix_get(vector, 0, 0, &entry);
    fixed_point_t norm = euclidean_norm(vector);
    for (int8_t i = 0; i < vector->rows; ++i) {
        matrix_get(vector, i, 0, &entry);
        matrix_set(vector, i, 0, fixed_point_divide(entry, norm));
    }
}

void gradient_ascent(matrix_t* matrix, matrix_t* vector) {
    matrix_t gradient;
    fixed_point_t vector_1_entry = fixed_point_init(0);
    fixed_point_t vector_2_entry = fixed_point_init(0);

    fixed_point_t previous_norm = fixed_point_init(0);
    for (int16_t i = 0; i < 1000; ++i) {
        matrix_multiply(matrix, vector, &gradient);
        fixed_point_t norm_AtAv = euclidean_norm(&gradient);

        if (fixed_point_absolute_value(
                fixed_point_subtract(norm_AtAv, previous_norm)) <
            fixed_point_divide(fixed_point_init(1), fixed_point_init(100))) {
            break;
        }
        previous_norm = norm_AtAv;

        for (int8_t j = 0; j < vector->rows; ++j) {
            matrix_get(vector, j, 0, &vector_1_entry);
            matrix_get(&gradient, j, 0, &vector_2_entry);  // gradient
            matrix_set(
                vector, j, 0,
                fixed_point_add(
                    vector_1_entry,
                    fixed_point_divide(vector_2_entry, fixed_point_init(128))));
        }

        normalize(vector);
    }
}

bool svd_calculate_v(matrix_t* result) {
    matrix_t singular_vector, matrix_column, outer_product_matrix,
        intermediate_matrix, second_vector_computing_matrix;
    matrix_init(&singular_vector, AtA.cols, 1);
    matrix_init(&matrix_column, AtA.cols, 1);
    matrix_init(&outer_product_matrix, AtA.rows, AtA.cols);
    matrix_init(&intermediate_matrix, AtA.rows, AtA.cols);
    matrix_init(&second_vector_computing_matrix, AtA.rows, AtA.cols);

    for (int8_t i = 0; i < AtA.cols; ++i) {
        if (i == 0) {
            for (int8_t j = 0; j < singular_vector.rows; ++j) {
                matrix_set(&singular_vector, j, 0, fixed_point_init(1));
            }
            gradient_ascent(&AtA, &singular_vector);
        } else if (i == 1) {
            matrix_set(&singular_vector, 0, 0, fixed_point_init(1));
            matrix_set(&singular_vector, 1, 0, -fixed_point_init(2));
            matrix_set(&singular_vector, 2, 0, fixed_point_init(1));

            matrix_get_column(result, 0, &matrix_column);
            identity_subtract_outer_product(&matrix_column,
                                            &outer_product_matrix);
            matrix_multiply(&outer_product_matrix, &AtA, &intermediate_matrix);
            matrix_multiply(&intermediate_matrix, &outer_product_matrix,
                            &second_vector_computing_matrix);

            gradient_ascent(&second_vector_computing_matrix, &singular_vector);
        } else {
            cross_product(result, &singular_vector);
        }

        normalize(&singular_vector);

        fixed_point_t entry = fixed_point_init(0);
        for (int8_t j = 0; j < result->rows; ++j) {
            matrix_get(&singular_vector, j, 0, &entry);
            matrix_set(result, j, i, entry);
        }
    }
    return true;
}
