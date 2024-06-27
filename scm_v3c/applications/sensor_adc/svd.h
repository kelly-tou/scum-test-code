
// ONLY FOR 2X2 MATRICES

#ifndef __SVD_H
#define __SVD_H

#include <stdbool.h>
#include <stdint.h>

#include "matrix.h"

bool svd_init(matrix_t* mtx);

// right singular vectors matrix
void svd_calculate_v(matrix_t* result);

#endif // __SVD_H
