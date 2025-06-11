#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "fixed_point.h"
#include "matrix.h"

matrix_t AtA;
fixed_point_t singular_vector_buffer[3], matrix_column_buffer[3],
    outer_product_matrix_buffer[9], intermediate_matrix_buffer[9],
    second_vector_computing_matrix_buffer[9], gradient_buffer[3];

bool svd_init(const matrix_t* matrix) {
    AtA = *matrix;
    return true;
}

void matrix_get_column(const matrix_t* matrix, const int8_t col,
                       matrix_t* result) {
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
        matrix_multiply(matrix, vector, &gradient, gradient_buffer);
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
    for (int8_t j = 0; j < result->rows; ++j) {
        matrix_set(result, j, 0, fixed_point_init(1));
    }
    gradient_ascent(&AtA, result);

    normalize(&result);

    return true;
}