#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "adc.h"
#include "fixed_point.h"
#include "matrix.h"
#include "memory_map.h"
#include "optical.h"
#include "scm3c_hw_interface.h"
#include "svd.h"

// Number of for loop cycles between ADC reads.
// 700000 for loop cycles roughly correspond to 1 second.
#define NUM_CYCLES_BETWEEN_ADC_READS 700000

matrix_t g_matrix, v;
fixed_point_short_t a, b, c, d;

// ADC configuration.
static const adc_config_t g_adc_config = {
    .reset_source = ADC_RESET_SOURCE_FSM,
    .convert_source = ADC_CONVERT_SOURCE_FSM,
    .pga_amplify_source = ADC_PGA_AMPLIFY_SOURCE_FSM,
    .pga_gain = 0,
    .settling_time = 0,
    .bandgap_reference_tuning_code = 1,
    .const_gm_tuning_code = 0xFF,
    .vbat_div_4_enabled = false,
    .ldo_enabled = true,
    .input_mux_select = ADC_INPUT_MUX_SELECT_EXTERNAL_SIGNAL,
    .pga_bypass = true,
};

int main(void) {
    initialize_mote();

    // Configure the ADC.
    printf("Configuring the ADC.\n");
    adc_config(&g_adc_config);
    adc_enable_interrupt();

    analog_scan_chain_write();
    analog_scan_chain_load();

    crc_check();
    perform_calibration();

    while (true) {
        // printf("Reading the ADC output.\n");
        // uint16_t adc_output = adc_read_output();
        // printf("ADC output: %u\n", adc_output);

        matrix_init(&g_matrix, 2, 2);
        matrix_init(&v, 2, 2);
        matrix_set(&g_matrix, 0, 0, fixed_point_init(4));   // a
        matrix_set(&g_matrix, 0, 1, fixed_point_init(0));   // b
        matrix_set(&g_matrix, 1, 0, fixed_point_init(3));   // c
        matrix_set(&g_matrix, 1, 1, fixed_point_init(-5));  // d

        // matrix_set(&g_matrix, 0, 0, fixed_point_init(1)); // a
        // matrix_set(&g_matrix, 0, 1, fixed_point_init(0)); // b
        // matrix_set(&g_matrix, 1, 0, fixed_point_init(0)); // c
        // matrix_set(&g_matrix, 1, 1, fixed_point_init(1)); // d

        svd_init(&g_matrix);
        svd_calculate_v(&v);

        matrix_get(&v, 0, 0, &a);
        matrix_get(&v, 0, 1, &b);
        matrix_get(&v, 1, 0, &c);
        matrix_get(&v, 1, 1, &d);

        // matrix_get(&g_matrix, 0, 0, &a);
        // matrix_get(&g_matrix, 0, 1, &b);
        // matrix_get(&g_matrix, 1, 0, &c);
        // matrix_get(&g_matrix, 1, 1, &d);

        // a = fixed_point_init(4);
        // b = fixed_point_init(0);
        // c = fixed_point_init(3);
        // d = fixed_point_init(-5);

        printf("a = %d, b = %d, c = %d, d = %d\n", a, b, c, d);

        // Wait for the next ADC read.
        for (size_t i = 0; i < NUM_CYCLES_BETWEEN_ADC_READS; ++i) {}
    }
}
