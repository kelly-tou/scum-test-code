#ifndef __LINEAR_REGRESSION_H
#define __LINEAR_REGRESSION_H

#include <stdbool.h>
#include <stdint.h>

#include "fixed_point.h"

bool linear_regression_init(const fixed_point_t signal_sampling_period,
                            const uint16_t* data, const size_t length);

bool linear_regression_get_time_constant(fixed_point_t* time_constant);

#endif
