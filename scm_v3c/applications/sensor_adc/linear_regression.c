#include "linear_regression.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "fixed_point.h"

#define MAX_NUM_SAMPLES 5000

fixed_point_t A[MAX_NUM_SAMPLES], B[MAX_NUM_SAMPLES];
fixed_point_t A_magnitude, A_B_dot_product;

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

bool linear_regression_init(const fixed_point_t signal_sampling_period,
                            const uint16_t* data, const size_t length) {
    fixed_point_t a_entry = fixed_point_init(0);
    fixed_point_t b_entry = fixed_point_init(0);
    fixed_point_t minimum_value = fixed_point_init(data[length - 1]);
    fixed_point_t amplitude = fixed_point_init(data[0]);
    A_magnitude = fixed_point_init(0);
    A_B_dot_product = fixed_point_init(0);

    for (int i = 0; i < length; ++i) {
        a_entry = fixed_point_divide(i, signal_sampling_period);
        A[i] = a_entry;
        A_magnitude = fixed_point_add(A_magnitude,
                                      fixed_point_multiply(a_entry, a_entry));
        if (!ln(fixed_point_divide(fixed_point_init(data[i]) - minimum_value,
                                   amplitude),
                b_entry)) {
            return false;
        }
        B[i] = -b_entry;
        A_B_dot_product = fixed_point_add(
            A_B_dot_product, fixed_point_multiply(a_entry, b_entry));
    }

    return true;
}

bool linear_regression_get_time_constant(fixed_point_t* time_constant) {
    *time_constant = fixed_point_divide(A_magnitude, A_B_dot_product);
    return true;
}
