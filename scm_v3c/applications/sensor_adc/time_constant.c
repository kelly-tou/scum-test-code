#include "time_constant.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "adc_msb.h"
#include "fixed_point.h"
#include "linear_regression.h"
// #include "matrix_pencil_method.h"

// Maximum number of ADC samples.
#define TIME_CONSTANT_MAX_NUM_ADC_SAMPLES 5000

// Number of samples to average at the end to find the minimum ADC sample.
#define TIME_CONSTANT_NUM_AVERAGES_FOR_MIN_ADC_SAMPLE 100

// Empirical ADC standard deviation in LSBs.
#define TIME_CONSTANT_ADC_STDDEV 5

// Time constant analysis state enum.
typedef enum {
    TIME_CONSTANT_STATE_INVALID = -1,
    TIME_CONSTANT_STATE_IDLE = 0,
    TIME_CONSTANT_STATE_RECEIVING = 1,
    TIME_CONSTANT_STATE_HAS_SUFFICIENT_SAMPLES = 2,
    TIME_CONSTANT_STATE_DONE = 3,
} time_constant_state_t;

// Time constant measurement frequency in Hz.
static uint32_t g_time_constant_sampling_frequency = 0;

// Time constant state.
static time_constant_state_t g_time_constant_state =
    TIME_CONSTANT_STATE_INVALID;

// Number of ADC samples in the buffer.
static size_t g_time_constant_num_adc_samples = 0;

// Buffer for the ADC samples.
static uint16_t g_time_constant_adc_samples[TIME_CONSTANT_MAX_NUM_ADC_SAMPLES];

// Minimum ADC sample in the buffer.
static uint16_t g_time_constant_min_adc_sample = ADC_MAX_THEORETICAL_ADC_SAMPLE;

// Maximum ADC sample in the buffer.
static uint16_t g_time_constant_max_adc_sample = 0;

// Number of samples since the minimum ADC sample was updated.
static size_t g_time_constant_num_samples_since_minimum = 0;

// NEW
static size_t g_time_constant_max_num_samples_since_minimum = 0;

// Add an ADC sample.
static inline void time_constant_add_to_buffer(const uint16_t sample) {
    g_time_constant_adc_samples[g_time_constant_num_adc_samples] = sample;
    ++g_time_constant_num_adc_samples;
}

void time_constant_init(const uint32_t time_constant_sampling_period_ms) {
    g_time_constant_sampling_frequency =
        1000 / time_constant_sampling_period_ms;
    g_time_constant_num_adc_samples = 0;
    g_time_constant_min_adc_sample = ADC_MAX_THEORETICAL_ADC_SAMPLE;
    g_time_constant_num_samples_since_minimum = 0;
    g_time_constant_max_num_samples_since_minimum = 0;
    g_time_constant_state = TIME_CONSTANT_STATE_IDLE;
}

bool time_constant_add_sample(const uint16_t adc_sample) {
    switch (g_time_constant_state) {
        case TIME_CONSTANT_STATE_IDLE: {
            const uint16_t previous_adc_sample =
                g_time_constant_adc_samples[g_time_constant_num_adc_samples];
            if (previous_adc_sample == ADC_MAX_ADC_SAMPLE &&
                adc_sample < ADC_MAX_ADC_SAMPLE) {
                adc_msb_init(adc_sample + ADC_MAX_THEORETICAL_ADC_SAMPLE -
                             ADC_MAX_ADC_SAMPLE);
                const uint16_t disambiguated_adc_sample =
                    adc_msb_get_disambiguated_sample();
                time_constant_add_to_buffer(disambiguated_adc_sample);
                g_time_constant_state = TIME_CONSTANT_STATE_RECEIVING;
            } else {
                g_time_constant_adc_samples[g_time_constant_num_adc_samples] =
                    adc_sample;
            }
            break;
        }
        case TIME_CONSTANT_STATE_RECEIVING: {
            adc_msb_disambiguate(adc_sample);
            const uint16_t disambiguated_adc_sample =
                adc_msb_get_disambiguated_sample();
            time_constant_add_to_buffer(disambiguated_adc_sample);

            // Update the minimum ADC sample.
            if (disambiguated_adc_sample <= 128 &&
                g_time_constant_min_adc_sample ==
                    ADC_MAX_THEORETICAL_ADC_SAMPLE) {
                g_time_constant_min_adc_sample = disambiguated_adc_sample;
                g_time_constant_num_samples_since_minimum = 0;
                g_time_constant_max_num_samples_since_minimum =
                    4 * g_time_constant_num_adc_samples / 5;
            } else if (g_time_constant_min_adc_sample !=
                       ADC_MAX_THEORETICAL_ADC_SAMPLE) {
                ++g_time_constant_num_samples_since_minimum;

                // If the minimum ADC sample has not been updated for a while,
                // the exponential might have finished decaying.
                if (g_time_constant_num_samples_since_minimum >
                    g_time_constant_max_num_samples_since_minimum) {
                    g_time_constant_state =
                        TIME_CONSTANT_STATE_HAS_SUFFICIENT_SAMPLES;
                }
            }
            break;
        }
        case TIME_CONSTANT_STATE_INVALID:
        case TIME_CONSTANT_STATE_HAS_SUFFICIENT_SAMPLES:
        case TIME_CONSTANT_STATE_DONE:
        default: {
            return false;
        }
    }
    return true;
}

bool time_constant_has_sufficient_samples(void) {
    return g_time_constant_state == TIME_CONSTANT_STATE_HAS_SUFFICIENT_SAMPLES;
}

fixed_point_t time_constant_estimate(void) {
    // matrix_pencil_method_init(
    //     fixed_point_init(g_time_constant_sampling_frequency),
    //     g_time_constant_adc_samples, g_time_constant_num_adc_samples);
    // fixed_point_t time_constant = fixed_point_init(0);
    // matrix_pencil_method_get_time_constant(&time_constant);

    linear_regression_init(fixed_point_init(g_time_constant_sampling_frequency),
                           g_time_constant_adc_samples,
                           g_time_constant_num_adc_samples);
    fixed_point_t time_constant = fixed_point_init(0);
    linear_regression_get_time_constant(&time_constant);

    return time_constant;
}
