
// // ONLY FOR MX3 MATRICES

// #include "svd_2.h"

// #include <stdbool.h>
// #include <stddef.h>
// #include <stdint.h>

// #include "fixed_point.h"
// #include "matrix.h"

// matrix_t svd_computing_matrix;
// matrix_t transposed_matrix;

// void matrix_print(const matrix_t* matrix) {
//     fixed_point_t entry = fixed_point_init(0);
//     for (size_t i = 0; i < matrix->rows; ++i) {
//         for (size_t j = 0; j < matrix->cols; ++j) {
//             matrix_get(matrix, i, j, &entry);
//             printf("%lld  ", entry);
//         }
//         printf("\n");
//     }
//     printf("\n");
// }

// // Using the Newton-Raphson Method, find the first root of a cubic polynomial
// of the form: x^3 + bx^2 + cx + d = 0.
// // Return whether finding the first root was successful.
// static bool find_cubic_root(const fixed_point_t b, const fixed_point_t c,
// const fixed_point_t d,
//                             const fixed_point_t initial_guess, fixed_point_t*
//                             root) {
//     fixed_point_t x = initial_guess;
//     for (int16_t i = 0; i < 1000; ++i) {
//         fixed_point_t fx =
//         fixed_point_add(fixed_point_add(fixed_point_cube(x),
//                                                             fixed_point_multiply(b,
//                                                             fixed_point_square(x))),
//                                             fixed_point_add(fixed_point_multiply(c,
//                                             x),
//                                                             d));
//         fixed_point_t fpx =
//         fixed_point_add(fixed_point_add(fixed_point_multiply(fixed_point_init(3),
//                                                                                 fixed_point_square(x)),
//                                                             fixed_point_multiply(fixed_point_init(2),
//                                                                                 fixed_point_multiply(x, b))),
//                                             c);
//         if (fpx == fixed_point_init(0)){
//             return false;
//         }
//         fixed_point_t x_new = fixed_point_subtract(x, fixed_point_divide(fx,
//         fpx)); if
//         (fixed_point_multiply(fixed_point_absolute_value(fixed_point_subtract(x_new,
//         x)), fixed_point_init(10000000)) < fixed_point_init(1)) {
//             *root = x_new;
//             return true;
//         }
//         x = x_new;
//     }
//     *root = x;
//     return true;
// }

// // Perform synthetic division to divide a cubic polynomial by a linear factor
// for cubic polynomials of the form: x^3 + bx^2 + cx + d = 0.
// // Return whether synthetic division was successful.
// static bool synthetic_division(const fixed_point_t b, const fixed_point_t c,
// fixed_point_t* b_new,
//                                fixed_point_t* c_new, const fixed_point_t
//                                root) {
//     *b_new = fixed_point_add(b, root);
//     *c_new = fixed_point_add(c, fixed_point_multiply(root, fixed_point_add(b,
//     root))); return true;
// }

// // Using the Jenkins Traub Algorithm, find the roots of a cubic polynomial of
// the form: x^3 + bx^2 + cx + d = 0.
// // Return whether finding the roots was successful.
// static bool jenkins_traub_cubic(const fixed_point_t b, const fixed_point_t c,
// const fixed_point_t d,
//                                 fixed_point_t* lambda1, fixed_point_t*
//                                 lambda2, fixed_point_t* lambda3) {
//     if (find_cubic_root(b, c, d, fixed_point_init(0), lambda1)) {
//         fixed_point_t b_new = fixed_point_init(0);
//         fixed_point_t c_new = fixed_point_init(0);
//         if (synthetic_division(b, c, &b_new, &c_new, *lambda1)) {
//             fixed_point_t discriminant =
//             fixed_point_subtract(fixed_point_square(b_new),
//             fixed_point_multiply(fixed_point_init(4), c_new)); if
//             (discriminant >= fixed_point_init(0)) {
//                 *lambda2 =
//                 fixed_point_divide(fixed_point_subtract(fixed_point_square_root(discriminant),
//                 b_new), fixed_point_init(2)); *lambda3 =
//                 fixed_point_divide(fixed_point_subtract(-fixed_point_square_root(discriminant),
//                 b_new), fixed_point_init(2)); return true;
//             }
//             return false;
//         }
//         return false;
//     }
//     return false;
// }

// // In the given matrix, replace the row at row_index2 with (constant * the
// row at row_index1) + the row at row_index2.
// // Return whether the replacement was successful.
// static bool add_multiples_of_rows(matrix_t* matrix, const size_t row_index1,
// const size_t row_index2, const size_t leading_entry_index,
//                            const fixed_point_t constant) {
//     fixed_point_t element_row1 = fixed_point_init(0);
//     fixed_point_t element_row2 = fixed_point_init(0);
//     for (size_t i = leading_entry_index; i < matrix->cols; ++i) {
//         matrix_get(matrix, row_index1, i, &element_row1);
//         matrix_get(matrix, row_index2, i, &element_row2);
//         matrix_set(matrix, row_index2, i, fixed_point_add(element_row2,
//         fixed_point_multiply(element_row1, constant)));
//     }
//     return true;
// }

// // Find and return the column index of the leading entry of the specified row
// of the given matrix. static size_t get_leading_entry_index(const matrix_t*
// matrix, const size_t row) {
//     fixed_point_t element = fixed_point_init(0);
//     for (size_t i = 0; i < matrix->cols; ++i) {
//         matrix_get(matrix, row, i, &element);
//         if (fixed_point_integer(element) != 0) {
//             return i;
//         }
//     }
//     return -1;
// }

// // Sort and modify a list of 3 elements in descending order.
// static bool sort_list_descending_order_and_determine_duplicates(fixed_point_t
// arr[3]) {
//     fixed_point_t temp = fixed_point_init(0);
//     // Compare and swap first two elements if needed
//     if (arr[0] < arr[1]) {
//         temp = arr[0];
//         arr[0] = arr[1];
//         arr[1] = temp;
//     }
//     // Compare and swap first and third elements if needed
//     if (arr[0] < arr[2]) {
//         temp = arr[0];
//         arr[0] = arr[2];
//         arr[2] = temp;
//     }
//     // Compare and swap second and third elements if needed
//     if (arr[1] < arr[2]) {
//         temp = arr[1];
//         arr[1] = arr[2];
//         arr[2] = temp;
//     }
//     if (arr[0] == arr[1] || arr[0] == arr[2] || arr[1] == arr[2]) {
//         return true;
//     }
//     return false;
// }

// bool matrix_eigenvalues(const matrix_t* matrix, fixed_point_t* lambda1,
// fixed_point_t* lambda2, fixed_point_t* lambda3) {
//     if (matrix->rows == 3 && matrix->cols == 3) {
//         fixed_point_t a = fixed_point_init(0);
//         fixed_point_t b = fixed_point_init(0);
//         fixed_point_t c = fixed_point_init(0);
//         fixed_point_t d = fixed_point_init(0);
//         fixed_point_t e = fixed_point_init(0);
//         fixed_point_t f = fixed_point_init(0);
//         fixed_point_t g = fixed_point_init(0);
//         fixed_point_t h = fixed_point_init(0);
//         fixed_point_t k = fixed_point_init(0);

//         matrix_get(matrix, 0, 0, &a);
//         matrix_get(matrix, 0, 1, &b);
//         matrix_get(matrix, 0, 2, &c);
//         matrix_get(matrix, 1, 0, &d);
//         matrix_get(matrix, 1, 1, &e);
//         matrix_get(matrix, 1, 2, &f);
//         matrix_get(matrix, 2, 0, &g);
//         matrix_get(matrix, 2, 1, &h);
//         matrix_get(matrix, 2, 2, &k);

//         fixed_point_t lambda_squared_coefficient = -fixed_point_add(a,
//         fixed_point_add(e, k)); fixed_point_t lambda_coefficient =
//         fixed_point_subtract(fixed_point_multiply(lambda_squared_coefficient,
//         fixed_point_divide(lambda_squared_coefficient, fixed_point_init(2))),
//                                                                 fixed_point_add(fixed_point_add(fixed_point_add(fixed_point_add(fixed_point_add(fixed_point_multiply(a,
//                                                                 fixed_point_divide(a,
//                                                                 fixed_point_init(2))),
//                                                                                                                                                 fixed_point_multiply(e, fixed_point_divide(e, fixed_point_init(2)))),
//                                                                                                                                 fixed_point_multiply(k, fixed_point_divide(k, fixed_point_init(2)))),
//                                                                                                                 fixed_point_multiply(b, d)),
//                                                                                                 fixed_point_multiply(c, g)),
//                                                                                 fixed_point_multiply(f, h)));
//         fixed_point_t constant =
//         -fixed_point_add(fixed_point_add(fixed_point_multiply(a,
//         fixed_point_subtract(fixed_point_multiply(e, k),
//                                                                                                                 fixed_point_multiply(f, h))),
//                                                                   fixed_point_multiply(b,
//                                                                   fixed_point_subtract(fixed_point_multiply(f,
//                                                                   g),
//                                                                                                                 fixed_point_multiply(d, k)))),
//                                                   fixed_point_multiply(c,
//                                                   fixed_point_subtract(fixed_point_multiply(d,
//                                                   h),
//                                                                                                 fixed_point_multiply(e, g))));
//         if (jenkins_traub_cubic(lambda_squared_coefficient,
//         lambda_coefficient, constant, lambda1, lambda2, lambda3)) {
//             return true;
//         }
//         return false;
//     }
//     return false;
// }

// bool matrix_eigenvectors(const matrix_t* matrix, const fixed_point_t lambda1,
// const fixed_point_t lambda2, const fixed_point_t lambda3, matrix_t* result) {
//     if (result->rows != 3 || result->cols != 3 || matrix-> rows != 3 ||
//     matrix->cols != 3){
//         return false;
//     }
//     fixed_point_t eigenvalues[3] = {lambda1, lambda2, lambda3};
//     bool repeated_eigenvalues =
//     sort_list_descending_order_and_determine_duplicates(eigenvalues);

//     matrix_t temp_matrix;
//     fixed_point_t diagonal_element = fixed_point_init(0);
//     fixed_point_t leading_entry = fixed_point_init(0);
//     fixed_point_t column_entry = fixed_point_init(0);
//     fixed_point_t temp_matrix_entry = fixed_point_init(0);
//     size_t leading_entry_index = 0;

//     for (size_t eigenvalue_index = 0; eigenvalue_index < matrix->rows;
//     ++eigenvalue_index) {
//         // printf("Eigenvalue: %lld\n", eigenvalues[eigenvalue_index]);
//         // subtract the eigenvalue from the diagonal entry
//         matrix_copy(matrix, &temp_matrix);
//         // printf("Start: \n");
//         matrix_print(&temp_matrix);
//         for (size_t diagonal_index = 0; diagonal_index < temp_matrix.rows;
//         ++diagonal_index) {
//             matrix_get(&temp_matrix, diagonal_index, diagonal_index,
//             &diagonal_element); matrix_set(&temp_matrix, diagonal_index,
//             diagonal_index, fixed_point_subtract(diagonal_element,
//             eigenvalues[eigenvalue_index]));
//         }

//         // matrix_print(&temp_matrix);

//         // row reduce
//         for (size_t current_row = 0; current_row < temp_matrix.rows - 1;
//         ++current_row) {
//             leading_entry_index = get_leading_entry_index(&temp_matrix,
//             current_row); if (leading_entry_index != -1) {
//                 matrix_get(&temp_matrix, current_row, leading_entry_index,
//                 &leading_entry); for (size_t other_row = 0; other_row <
//                 temp_matrix.rows; ++other_row) {
//                     if (other_row != current_row) {
//                         matrix_get(&temp_matrix, other_row,
//                         leading_entry_index, &column_entry);
//                         add_multiples_of_rows(&temp_matrix, current_row,
//                         other_row, leading_entry_index,
//                         -fixed_point_divide(column_entry, leading_entry));
//                     }
//                 }
//             }
//             // matrix_print(&temp_matrix);
//         }

//         matrix_set(&temp_matrix, 2, 0, fixed_point_init(0));
//         matrix_set(&temp_matrix, 2, 1, fixed_point_init(0));
//         matrix_set(&temp_matrix, 2, 2, fixed_point_init(0));

//         // matrix_print(&temp_matrix);

//         // determine the free variables
//         size_t free_variable_index = -1;
//         for (size_t row = 0; row < temp_matrix.rows; ++row) {
//             if (repeated_eigenvalues) {
//                 leading_entry_index = get_leading_entry_index(&temp_matrix,
//                 temp_matrix.rows - row - 1);
//             } else {
//                 leading_entry_index = get_leading_entry_index(&temp_matrix,
//                 row);
//             }
//             if (leading_entry_index == -1) {
//                 if (repeated_eigenvalues) {
//                     free_variable_index = temp_matrix.rows - row - 1;
//                 } else {
//                     free_variable_index = row;
//                 }
//                 break;
//             }
//         }

//         // find the right singular vector corresponding to this eigenvalue
//         fixed_point_t norm = fixed_point_init(1);
//         matrix_set(result, free_variable_index, eigenvalue_index, norm);
//         for (size_t row = 0; row < temp_matrix.rows; ++row) {
//             if (row != free_variable_index) {
//                 leading_entry_index = get_leading_entry_index(&temp_matrix,
//                 row); matrix_get(&temp_matrix, row, leading_entry_index,
//                 &leading_entry); matrix_get(&temp_matrix, row,
//                 free_variable_index, &column_entry); fixed_point_t entry =
//                 -fixed_point_divide(column_entry, leading_entry);
//                 matrix_set(result, row, eigenvalue_index, entry);
//                 norm = fixed_point_add(norm, fixed_point_square(entry));
//             }
//         }
//         norm = fixed_point_square_root(norm);

//         // printf("before normalize: \n");
//         // matrix_print(result);

//         // normalize the right singular vector
//         for (size_t vector_row_index = 0; vector_row_index < result->rows;
//         ++vector_row_index) {
//             matrix_get(result, vector_row_index, eigenvalue_index,
//             &column_entry); matrix_set(result, vector_row_index,
//             eigenvalue_index, fixed_point_divide(column_entry, norm));
//         }

//         // printf("after normalize: \n");
//         // matrix_print(result);
//     }
//     return true;
// }

// bool svd_init(const matrix_t* matrix) {
//     if (matrix->cols != 3) {
//         return false;
//     }
//     matrix_transpose(matrix, &transposed_matrix);
//     matrix_multiply(&transposed_matrix, matrix, &svd_computing_matrix);
//     return true;
// }

// bool svd_init_given_transpose_product(const matrix_t* matrix) {
//     svd_computing_matrix = *matrix;
//     // printf("data matrix: \n");
//     // matrix_print(&svd_computing_matrix);
//     return true;
// }

// bool svd_calculate_v(matrix_t* result) {
//     fixed_point_t lambda1 = fixed_point_init(0);
//     fixed_point_t lambda2 = fixed_point_init(0);
//     fixed_point_t lambda3 = fixed_point_init(0);
//     matrix_eigenvalues(&svd_computing_matrix, &lambda1, &lambda2, &lambda3);
//     matrix_eigenvectors(&svd_computing_matrix, lambda1, lambda2, lambda3,
//     result); return true;
// }
