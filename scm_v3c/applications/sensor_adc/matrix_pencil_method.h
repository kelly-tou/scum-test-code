#ifndef __MATRIX_PENCIL_METHOD_H
#define __MATRIX_PENCIL_METHOD_H

#include <stdbool.h>
#include <stdint.h>

#include "fixed_point.h"
#include "matrix.h"
#include "svd_3_v2.h"

bool matrix_pencil_method_init(const fixed_point_t signal_sampling_period,
                               int16_t data[], const size_t length);

bool matrix_pencil_method_get_time_constant(fixed_point_t* time_constant);

#endif
