
// ONLY FOR 2X2 MATRICES

#ifndef __SVD_H
#define __SVD_H

#include <stdbool.h>
#include <stdint.h>

#include "matrix.h"

// A = U(Lambda)(V^T)
// Initialize an SVD computation for the given matrix.
// Return whether the matrix was successfully intialized.
bool svd_init(matrix_t* matrix);

// Compute the V matrix and write the result to another matrix.
// Return whether the computation was successful.
void svd_calculate_v(matrix_t* result);

#endif  // __SVD_H
