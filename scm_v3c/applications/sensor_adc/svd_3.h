// calculate the v matrix of the SVD for 3x3 matrices only

#ifndef __SVD_3_H
#define __SVD_3_H

#include <stdbool.h>
#include <stdint.h>

#include "fixed_point.h"
#include "matrix.h"

// matrix = (A^T)A
// Initialize an SVD computation for the given matrix.
// Return whether the matrix was successfully intialized.
bool svd_init(const matrix_t* matrix);

// Compute the V matrix and write the result to another matrix.
// Return whether the computation was successful.
bool svd_calculate_v(matrix_t* result);

#endif  // __SVD_3_H
