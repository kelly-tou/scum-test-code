#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "adc.h"
#include "fixed_point.h"
#include "matrix.h"
#include "matrix_pencil_method.h"
#include "memory_map.h"
#include "optical.h"
#include "scm3c_hw_interface.h"

// Number of for loop cycles between ADC reads.
// 700000 for loop cycles roughly correspond to 1 second.
#define NUM_CYCLES_BETWEEN_ADC_READS 700000

int16_t adc_samples[] = {
    1023, 950, 868, 798, 736, 672, 623, 587, 548, 511, 497, 468, 430, 407, 383,
    359,  342, 319, 310, 287, 274, 255, 255, 242, 229, 214, 207, 202, 190, 191,
    190,  181, 171, 166, 161, 151, 156, 145, 143, 143, 142, 141, 127, 143, 135,
    127,  127, 127, 127, 127, 127, 127, 126, 124, 127, 126, 127, 127, 127, 127,
    120,  120, 115, 118, 120, 111, 120, 118, 118, 117, 114, 112, 111, 115, 119,
    119,  119, 120, 111, 120, 111, 118, 114, 103, 119, 111, 118, 115, 112, 110,
    112,  111, 122, 111, 115, 111, 119, 119, 111, 119, 118, 114, 111, 111, 119,
    114,  111, 118, 111, 119, 116, 111, 119, 111, 115, 114, 107, 118, 111, 114,
    103,  112, 110, 111, 111, 118, 111, 118, 110, 111, 118, 111, 119, 111, 111,
    118,  111, 115, 111, 111, 118, 111, 116, 107, 112, 103, 111, 119, 115, 111,
    111,  120, 111, 119, 114, 111, 114, 110, 110, 111, 119, 111, 115, 111, 114,
    112,  106, 111, 111, 112, 107, 111, 120, 110, 111, 119, 112, 111, 118, 111,
    119,  111, 111, 127, 111, 118, 111, 111, 119, 113, 111, 111, 119, 118, 111,
    119,  113, 111, 115, 109, 114, 106, 111, 119, 112, 103, 111, 119, 114, 111,
    113,  108, 111, 117, 109, 114, 110, 111, 120, 107, 111, 115, 111, 115, 111,
    119,  112, 110, 126, 111, 118, 111, 111, 119, 116, 111, 112, 103, 115, 111,
    120,  110, 111, 119, 118, 111, 111, 111, 122, 111, 111, 111, 111, 119, 114,
    111,  112, 108, 111, 111, 116, 111, 111, 114, 103, 111, 111, 111, 115, 111,
    111,  116, 111, 119, 118, 111, 113, 108, 112, 103, 111, 119, 111, 114, 111,
    112,  106, 103, 115, 110, 111, 118, 111, 119, 111, 111, 119, 111, 118, 111,
    111};

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

    printf("start\r\n\r\n");

    // while (true) {
    // printf("Reading the ADC output.\n");
    // uint16_t adc_output = adc_read_output();
    // printf("ADC output: %u\n", adc_output);

    matrix_pencil_method_init(fixed_point_init(100), adc_samples, 301);
    fixed_point_t time_constant = fixed_point_init(0);

    matrix_pencil_method_get_time_constant(&time_constant);
    printf("Time Constant: %lld\n", time_constant);

    printf("stop\r\n\r\n");

    // Wait for the next ADC read.
    //     for (size_t i = 0; i < NUM_CYCLES_BETWEEN_ADC_READS; ++i) {}
    // }
}

/*

int16_t adc_samples[] = {1016,  992,  995,  982,  976,  971,  956,  967,  935,
934,  937, 910,  934,  903,  893,  903,  878,  887,  864,  850,  856,  846, 856,
828,  826,  830,  807,  817,  790,  786,  791,  766,  786, 754,  751,  767, 743,
749,  728,  714,  728,  700,  711,  702, 694,  703,  679,  680,  670,  660, 672,
647,  656,  643,  632, 646,  627,  631,  626,  613,  624,  607,  611,  607, 599,
609, 583,  590,  590,  575,  587,  570,  569,  571,  560,  567,  558, 551,  555,
537,  559,  542,  536,  537,  524,  535,  511,  531, 527,  511,  537,  505, 511,
511,  508,  511,  503,  503,  511, 490,  497,  479,  480,  482,  471,  483, 464,
455,  479,  453, 463,  447,  447,  466,  438,  447,  440,  431,  447,  430, 431,
                        430,  418,  431,  415,  423,  415,  407,  423,  395,
403,  399, 394,  407,  386,  383,  396,  374,  399,  374,  381,  383,  372, 383,
370,  367,  371,  355,  373,  351,  357,  359,  343,  367, 344,  339,  351, 335,
352,  326,  328,  335,  327,  343,  326, 319,  339,  311,  335,  319,  319, 335,
309,  319,  310,  303, 319,  303,  319,  303,  300,  308,  287,  303,  292, 282,
312, 277,  295,  287,  279,  303,  271,  287,  282,  263,  288,  255, 286,  271,
264,  276,  255,  279,  255,  255,  279,  255,  255, 255,  255,  279,  255, 255,
255,  255,  266,  240,  243,  255, 248,  255,  247,  246,  255,  244,  255, 240,
236,  242,  223, 254,  232,  226,  239,  227,  241,  223,  229,  238,  215, 239,
                        223,  224,  226,  203,  231,  215,  215,  231,  209,
226,  206, 207,  230,  203,  223,  207,  206,  223,  207,  222,  203,  199, 222,
191,  215,  206,  191,  223,  199,  201,  191,  191,  220, 191,  203,  191, 191,
212,  186,  191,  200,  179,  207,  183, 192,  190,  183,  207,  183,  191, 192,
169,  199,  175,  183, 191,  178,  191,  184,  175,  191,  175,  191,  191, 181,
187, 172,  186,  167,  170,  183,  171,  183,  168,  161,  182,  162, 175,  168,
161,  182,  161,  175,  166,  164,  174,  159,  175, 164,  158,  178,  152, 167,
159,  159,  175,  159,  174,  159, 159,  178,  147,  159,  159,  157,  175, 151,
159,  161,  143, 173,  151,  160,  151,  151,  169,  145,  152,  158,  149, 168,
                        139,  151,  155,  143,  170,  143,  151,  158,  143,
168,  135, 148,  158,  143,  167,  143,  143,  159,  143,  161,  134,  142, 159,
143,  159,  143,  143,  159,  143,  159,  143,  143,  159, 139,  151,  135, 135,
159,  135,  151,  139,  137,  151,  135, 151,  138,  134,  151,  127,  151, 136,
127,  159,  127,  151, 142,  127,  159,  127,  143,  143,  127,  159,  135, 142,
136, 119,  155,  127,  143,  143,  127,  159,  127,  143,  143,  127, 155,  127,
127,  143,  123,  151,  127,  135,  142,  127,  155, 127,  135,  143,  127, 152,
119,  127,  147,  126,  148,  126, 127,  150,  127,  151,  127,  127,  150, 127,
143,  124,  127, 143,  127,  144,  122,  127,  147,  119,  138,  120,  120, 143,
                        126,  143,  127,  127,  143,  126,  127,  127,  127,
150,  119, 127,  127,  127,  150,  122,  127,  135,  119,  143,  121,  127, 127,
126,  143,  124,  127,  136,  116,  139,  126,  127,  135, 123,  143,  119, 127,
143,  124,  143,  119,  127,  143,  121, 143,  127,  127,  143,  126,  138, 118,
127,  143,  119,  139, 119,  127,  143,  119,  135,  126,  119,  142,  117, 138,
113, 119,  142,  119,  127,  127,  127,  143,  120,  127,  127,  127, 143,  118,
127,  127,  123,  135,  118,  127,  127,  126,  143, 125,  127,  127,  127, 143,
122,  127,  127,  127,  143,  116, 126,  127,  125,  143,  120,  119,  127, 127,
151,  119,  123, 127,  126,  136,  111,  127,  139,  115,  142,  120,  120, 127,
                        120,  127,  125,  119,  135,  112,  127,  124,  119,
135,  119, 131,  115,  119,  138,  111,  135,  119,  122,  127,  120,  127, 127,
125,  135,  119,  127,  127,  127,  142,  114,  127,  127, 119,  139,  117, 127,
127,  120,  127,  120,  126,  127,  125, 138,  111,  127,  127,  122,  135, 116,
119,  127,  122,  135, 115,  124,  127,  123,  139,  111,  127,  127,  119, 143,
119, 126,  127,  122,  131,  111,  126,  127,  124,  135,  120,  119, 135,  111,
127,  126,  120,  127,  119,  135,  119,  119,  127, 119,  127,  122,  119, 131,
111,  127,  124,  115,  127,  119, 127,  127,  124,  135,  111,  127,  127, 120,
127,  119,  127, 123,  119,  135,  116,  127,  127,  120,  127,  119,  127, 127,
                        119,  137,  111,  127,  127,  119,  135,  115,  122,
127,  122, 127,  119,  124,  127,  122,  135,  114,  119,  127,  122,  127, 119,
126,  127,  124,  135,  118,  119,  127,  122,  127,  122, 122,  127,  122, 127,
122,  119,  131,  111,  127,  123,  120, 127,  124,  127,  119,  120,  131, 111,
127,  127,  123,  135, 111,  127,  127,  119,  136,  110,  127,  127,  120, 127,
120, 125,  127,  119,  134,  111,  127,  127,  119,  136,  111,  127, 127,  120,
127,  122,  126,  127,  119,  135,  118,  122,  127, 120,  127,  119,  127, 127,
119,  135,  114,  119,  127,  122, 127,  119,  125,  127,  119,  135,  119, 127,
135,  111,  127, 124,  119,  127,  120,  127,  122,  119,  135,  118,  127, 123,
                        122,  127,  119,  127,  127,  126,  135,  113,  126,
119,  119, 135,  111,  127,  126,  124,  127,  119,  127,  127,  119,  127, 123,
127,  127,  127,  135,  111,  127,  127,  126,  127,  119, 127,  127,  122, 127,
119,  127,  127,  119,  127,  119,  127, 127,  127,  135,  118,  122,  127, 123,
135,  118,  120,  127, 119,  135,  115,  120,  127,  119,  127,  127,  126, 127,
122, 127,  119,  120,  127,  119,  127,  124,  120,  127,  119,  127, 127,  119,
127,  119,  127,  127,  127,  127,  124,  127,  127, 124,  130,  110,  127, 126,
119,  135,  111,  127,  127,  122, 127,  119,  127,  127,  120,  127,  123, 127,
127,  126,  127, 122,  119,  127,  126,  138,  111,  126,  127,  123,  138, 111,
                        126,  127,  119,  127,  122,  119,  127,  127,  135,
118,  119, 127,  120,  127,  123,  127,  127,  122,  127,  127,  126,  135, 112,
127,  127,  126,  127,  125,  127,  123,  127,  136,  110, 127,  126,  119, 127,
119,  127,  127,  125,  135,  111,  127, 125,  120,  127,  119,  127,  127, 123,
135,  119,  122,  126, 124,  135,  111,  127,  127,  119,  149,  111,  123, 127,
125, 135,  118,  119,  127,  125,  135,  118}; length: 1041

int16_t adc_samples[] = {814, 559, 399, 296, 223, 191, 159, 143, 127, 127, 127,
126, 118, 111, 119, 118, 111, 111, 116, 112, 103, 111, 114, 106, 103, 111, 111,
112, 103, 111, 111, 111, 111, 111, 111, 111, 111, 112, 103, 103, 111, 111, 111,
114, 103, 110, 111, 111, 111, 111, 111, 111, 111, 114, 108, 104, 110, 111, 111,
112, 104, 106, 104, 106, 110, 111, 113, 110, 111, 114, 106, 107, 112, 103, 112,
                        104, 104, 102, 111, 111, 111, 111, 111, 111, 112, 104,
103, 111, 107, 111, 114, 103, 111, 112, 103, 111, 111, 111, 111, 111, 110, 106,
106, 110, 106, 104, 103, 112, 103, 107, 111, 114, 103, 111, 112, 103, 111, 116,
111, 112, 103, 111, 111, 118, 111, 111, 111, 111, 112, 104, 103, 111, 113, 110,
111, 111, 111, 111, 114, 106, 104, 104, 103, 111, 111, 115, 111, 112, 104, 105,
                        110, 111, 118, 111, 111, 111, 111, 111, 111, 111, 111,
111, 111, 111, 114, 109, 111, 112, 103, 111, 111, 116, 110, 111, 111, 111, 113,
103, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 111, 112, 103,
111, 111, 107, 111, 113, 111, 111, 111, 118, 111, 126, 111, 111, 116, 103, 112,
106, 111, 111, 111, 111, 111, 111, 112, 103, 111, 119, 111, 111, 111, 116, 103,
                        111, 114, 106, 106};
    length: 229

int16_t adc_samples[] = {1014, 834, 684, 564, 491, 424, 368, 318, 287, 254, 226,
206, 191, 178, 159, 159, 146, 127, 131, 127, 127, 127, 127, 120, 112, 111, 119,
118, 112, 107, 111, 111, 111, 111, 111, 112, 106, 107, 103, 111, 116, 109, 116,
111, 112, 103, 111, 111, 111, 112, 106, 111, 111, 111, 111, 118, 107, 110, 112,
103, 111, 111, 111, 114, 108, 111, 111, 111, 114, 111, 111, 112, 106, 106, 110,
                        111, 111, 114, 103, 111, 112, 107, 111, 111, 112, 108,
105, 107, 110, 106, 111, 111, 111, 114, 106, 103, 111, 118, 110, 111, 111, 114,
103, 109, 111, 111, 112, 101, 107, 111, 108, 111, 111, 116, 108, 103, 107, 111,
111, 111, 112, 101, 111, 112, 103, 111, 111, 111, 111, 112, 103, 111, 115, 111,
111, 111, 112, 99, 103, 111, 109, 109, 111, 111, 111, 112, 103, 111, 112, 101,
                        105, 106, 111, 103, 110, 111, 111, 111, 112, 110, 111,
111, 111, 111, 110, 111, 111, 106, 110, 111, 112, 99, 103, 110, 111, 111, 112,
103, 111, 111, 112, 101, 110, 107, 103, 111, 112, 106, 106, 109, 107, 112, 104,
103, 109, 110, 111, 106, 103, 111, 111, 110, 111, 111, 109, 111, 108, 107, 111,
114, 104, 106, 103, 111, 112, 103, 109, 111, 110, 103, 104, 103, 103, 110, 103,
                        110, 108, 103, 111, 111, 112, 95, 111, 111, 111, 113,
103, 111, 116, 104, 103, 111, 117, 103, 109, 111, 111, 127, 111, 109, 107, 111,
111, 111, 111, 103, 111, 112, 99, 111, 103, 109, 110, 110, 111, 111, 111, 111,
111, 111, 111, 111, 111, 116, 103, 109, 111, 111, 111, 111, 103, 111, 111, 111,
111, 112, 99, 104, 104, 103, 111, 111, 116, 103, 111, 112, 102, 102, 107, 103,
                        111, 110, 111, 111, 111, 113, 111, 109, 111, 111, 116,
103, 110, 110, 110, 109, 107, 106, 110, 103, 108, 111, 111, 111, 111, 111, 111,
111, 110, 107, 103, 111, 111, 107, 103, 111, 111, 111, 111, 111, 114, 103, 111,
103, 108, 103, 112, 103, 110, 111, 112, 106, 103, 111, 108, 103, 111, 111, 114,
101, 107, 107, 107, 105, 108, 103, 111, 111, 111, 115, 107, 111, 111, 111, 107,
                        111, 111, 111, 105, 103, 111, 110, 111};
    length: 383

int16_t adc_samples[] = {1018, 1019, 1012, 1004, 994, 983, 998, 988, 980, 972,
971, 970, 964, 952, 959, 959, 956, 942, 938, 942, 934, 927, 926, 928, 914, 918,
901, 907, 899, 894, 895, 900, 882, 882, 878, 871, 882, 866, 864, 854, 862, 860,
848, 840, 840, 832, 827, 831, 832, 821, 821, 824, 810, 809, 808, 804, 807, 806,
794, 783, 791, 782, 779, 776, 762, 767, 774, 764, 760, 749, 747, 758, 750, 744,
                        732, 727, 730, 728, 720, 710, 716, 710, 714, 702, 710,
703, 703, 695, 688, 682, 684, 679, 678, 672, 670, 672, 663, 663, 666, 657, 654,
654, 654, 650, 643, 639, 646, 631, 639, 638, 631, 629, 632, 623, 623, 617, 622,
615, 613, 607, 608, 599, 605, 605, 602, 598, 592, 590, 587, 587, 585, 575, 583,
575, 584, 574, 575, 575, 575, 575, 568, 559, 570, 562, 559, 559, 556, 552, 553,
                        545, 544, 539, 559, 550, 541, 535, 527, 533, 530, 530,
530, 526, 511, 535, 523, 519, 520, 511, 527, 511, 511, 534, 511, 511, 511, 511,
511, 511, 511, 511, 511, 511, 511, 511, 498, 503, 503, 503, 502, 498, 495, 496,
488, 492, 483, 485, 469, 479, 478, 479, 478, 480, 468, 466, 458, 466, 456, 462,
466, 455, 458, 462, 463, 455, 455, 450, 447, 451, 447, 447, 453, 441, 441, 444,
                        431, 431, 432, 424, 430, 439, 423, 423, 431, 434, 422,
423, 418, 406, 421, 419, 415, 421, 416, 407, 407, 407, 405, 410, 399, 407, 408,
399, 401, 391, 395, 399, 396, 391, 395, 391, 391, 393, 378, 381, 387, 380, 382,
383, 383, 383, 383, 383, 379, 380, 378, 376, 368, 362, 368, 367, 367, 370, 359,
367, 368, 359, 359, 363, 359, 359, 358, 352, 350, 354, 352, 352, 343, 350, 346,
                        347, 343, 347, 350, 343, 342, 336, 328, 327, 334, 335,
338, 327, 334, 327, 333, 323, 327, 330, 324, 319, 328, 320, 319, 325, 319, 319,
319, 327, 315, 316, 318, 318, 319, 319, 319, 319, 327, 314, 310, 303, 316, 312,
302, 310, 303, 319, 311, 311, 301, 299, 301, 299, 298, 291, 295, 295, 291, 295,
295, 294, 287, 295, 288, 279, 287, 287, 295, 283, 285, 286, 282, 279, 279, 287,
                        279, 279, 278, 279, 271, 275, 271, 277, 271, 275, 267,
271, 272, 266, 266, 259, 263, 269, 263, 263, 255, 271, 263, 255, 270, 255, 255,
267, 255, 266, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 254,
250, 249, 254, 255, 255, 255, 255, 255, 245, 247, 254, 250, 247, 247, 250, 248,
240, 240, 238, 240, 238, 236, 239, 243, 239, 242, 231, 231, 240, 226, 232, 226,
                        230, 231, 239, 238, 237, 234, 227, 229, 230, 229, 231,
229, 226, 231, 223, 232, 223, 231, 225, 224, 215, 224, 223, 228, 219, 221, 223,
223, 223, 223, 231, 220, 222, 223, 219, 218, 214, 214, 215, 212, 211, 211, 207,
215, 207, 215, 215, 210, 201, 204, 208, 206, 210, 207, 207, 208, 206, 207, 207,
207, 207, 207, 202, 200, 191, 206, 200, 191, 207, 199, 191, 201, 197, 199, 199,
                        191, 195, 191, 203, 191, 194, 191, 199, 198, 194, 191,
199, 189, 197, 191, 199, 191, 191, 191, 199, 191, 191, 191, 191, 191, 192, 186,
191, 191, 191, 191, 191, 191, 191, 191, 191, 191, 191, 191, 188, 184, 179, 181,
183, 190, 186, 179, 190, 190, 186, 183, 183, 184, 176, 175, 183, 183, 183, 183,
175, 179, 183, 183, 176, 174, 173, 175, 183, 179, 180, 174, 175, 178, 175, 178,
                        171, 172, 167, 175, 179, 175, 175, 175, 167, 172, 167,
173, 171, 174, 168, 166, 174, 170, 165, 168, 166, 170, 161, 163, 167, 167, 166,
163, 167, 167, 167, 165, 160, 159, 168, 159, 167, 168, 162, 159, 167, 165, 166,
159, 159, 163, 162, 154, 162, 160, 157, 159, 180, 155, 159, 159, 159, 167, 161,
155, 159, 163, 159, 160, 156, 158, 158, 159, 159, 159, 162, 154, 151, 159, 159,
                        159, 159, 159, 161, 151, 159, 159, 159, 159, 159, 159,
158, 159, 159, 159, 156, 151, 158, 151, 159, 151, 153, 147, 150, 143, 150, 154,
151, 151, 158, 151, 151, 150, 147, 147, 139, 143, 156, 145, 142, 143, 151, 151,
145, 144, 143, 150, 143, 151, 147, 142, 143, 152, 143, 143, 143, 151, 150, 146,
143, 143, 141, 143, 143, 148, 144, 139, 140, 136, 143, 145, 143, 143, 143, 143,
                        144, 135, 143, 143, 142, 139, 141, 143, 143, 143, 148,
135, 143, 143, 143, 146, 135, 143, 139, 136, 135, 143, 143, 143, 143, 139, 127,
143, 142, 138, 136, 135, 135, 135, 135, 133, 136, 127, 139, 127, 143, 135, 131,
127, 143, 139, 135, 135, 138, 134, 127, 143, 135, 136, 127, 142, 127, 145, 127,
143, 138, 127, 143, 140, 127, 139, 127, 142, 127, 139, 127, 135, 127, 142, 127,
                        140, 127, 141, 127, 131, 127, 139, 127, 143, 127, 143,
127, 142, 127, 142, 127, 142, 127, 135, 127, 127, 139, 132, 127, 142, 131, 127,
135, 127, 127, 127, 127, 127, 139, 127, 127, 135, 127, 135, 127, 127, 127, 127,
127, 127, 127, 135, 127, 127, 135, 127, 127, 127, 135, 127, 127, 135, 127, 133,
127, 131, 120, 127, 134, 127, 127, 134, 127, 127, 127, 127, 127, 127, 127, 127,
                        127, 131, 127, 127, 127, 127, 127, 127, 143, 127, 127,
127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
127, 127, 127, 127, 127, 127, 127, 127, 135, 127, 127, 127, 127, 127, 127, 127,
127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
                        127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
123, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127,
127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 127, 126, 122, 120, 119, 127,
127, 124, 122, 120, 119, 127, 127, 125, 125, 127, 127, 126, 127, 127, 127, 122,
119, 127, 124, 119, 127, 123, 126, 127, 127, 127, 127, 123, 124, 123, 126, 124,
                        120, 120, 119, 122, 112, 119, 123, 122, 119, 126, 126,
126, 124, 124, 120, 113, 119, 124, 120, 117, 126, 127, 124, 126, 119, 123, 121,
118, 119, 126, 122, 119, 124, 117, 122, 124, 119, 119, 125, 120, 116, 117, 119,
127, 122, 111, 122, 119, 122, 119, 126, 119, 120, 117, 119, 124, 119, 125, 125,
126, 127, 122, 122, 114, 115, 119, 119, 117, 117, 122, 115, 117, 119, 119, 122,
                        118, 118, 115, 120, 118, 120, 111, 127, 126, 125, 119,
120, 118, 120, 115, 115, 118, 119, 126, 124, 114, 114, 115, 122, 111, 119, 124,
119, 119, 126, 122, 119, 123, 119, 124, 116, 117, 119, 117, 119, 118, 119, 119,
125, 122, 119, 127, 127, 123, 119, 122, 117, 115, 116, 111, 119, 119, 122, 117,
119, 119, 119, 120, 111, 125, 117, 119, 126, 119, 119, 115, 112, 110, 115, 120,
                        114, 111, 122, 120, 111, 119, 122, 115, 111, 118, 119,
120, 111, 120, 118, 116, 120, 120, 112, 111, 120, 119, 120, 111, 119, 119, 119,
118, 115, 111, 119, 119, 119, 119, 120, 111, 119, 116, 118, 118, 114, 110, 114,
114, 112, 111, 119, 119, 116, 111, 118, 116, 112, 111, 119, 124, 115, 111, 119,
118, 119, 117, 112, 108, 113, 111, 119, 111, 119, 119, 111, 119, 120, 111, 120,
                        111, 116, 115, 114, 112, 111, 119, 119, 119, 116, 111,
119, 119, 119, 120, 110, 114, 116, 118, 111, 127, 122, 114, 111, 120, 112, 111,
119, 111, 122, 113, 111, 115, 118, 112, 110, 111, 120, 111, 123, 119, 111, 124,
119, 118, 119, 112, 112, 110, 111, 118, 116, 111, 120, 111, 115, 111, 118, 111,
119, 115, 120, 111, 111, 119, 111, 119, 119, 114, 104, 111, 116, 111, 119, 114,
                        111, 119, 116, 111, 119, 119, 119, 119, 117, 114, 111,
119, 115, 111, 119, 118, 112, 111, 123, 118, 113, 111, 119, 111, 119, 120, 114,
115, 111, 117, 111, 119, 119, 118, 114, 111, 117, 112, 105, 103, 111, 115, 118,
113, 111, 122, 116, 110, 112, 103, 112, 107, 111, 111, 119, 111, 111, 119, 116,
99, 110, 115, 111, 112, 111, 111, 118, 111, 117, 119, 118, 115, 112, 110, 115,
                        111, 112, 107, 111, 118, 111, 120, 112, 103, 111, 119,
118, 111, 114, 111, 119, 111, 118, 118, 112, 106, 111, 116, 111, 122, 111, 115,
112, 111, 114, 111, 118, 111, 116, 112, 111, 118, 115, 111, 117, 115, 115, 111,
118, 114, 111, 118, 116, 111, 116, 111, 115, 111, 114, 116, 111, 111, 119, 115,
111, 111, 119, 119, 119, 115, 115, 111, 120, 111, 111, 114, 111, 119, 114, 111,
                        120, 111, 111, 115, 111, 111, 118, 111, 111, 115, 114,
111, 119, 113, 111, 115, 111, 120, 111, 111, 117, 118, 112, 111, 118, 112, 111,
111, 107, 111, 114, 111, 119, 116, 111, 119, 111, 112, 110, 116, 111, 118, 111,
111, 118, 111, 118, 111, 111, 119, 118, 114, 103, 111, 116, 110, 111, 114, 111,
111, 111}; length: 1560

int16_t adc_samples[] = {1023, 950, 868, 798, 736, 672, 623, 587, 548, 511, 497,
468, 430, 407, 383, 359, 342, 319, 310, 287, 274, 255, 255, 242, 229, 214, 207,
202, 190, 191, 190, 181, 171, 166, 161, 151, 156, 145, 143, 143, 142, 141, 127,
143, 135, 127, 127, 127, 127, 127, 127, 127, 126, 124, 127, 126, 127, 127, 127,
127, 120, 120, 115, 118, 120, 111, 120, 118, 118, 117, 114, 112, 111, 115, 119,
                        119, 119, 120, 111, 120, 111, 118, 114, 103, 119, 111,
118, 115, 112, 110, 112, 111, 122, 111, 115, 111, 119, 119, 111, 119, 118, 114,
111, 111, 119, 114, 111, 118, 111, 119, 116, 111, 119, 111, 115, 114, 107, 118,
111, 114, 103, 112, 110, 111, 111, 118, 111, 118, 110, 111, 118, 111, 119, 111,
111, 118, 111, 115, 111, 111, 118, 111, 116, 107, 112, 103, 111, 119, 115, 111,
                        111, 120, 111, 119, 114, 111, 114, 110, 110, 111, 119,
111, 115, 111, 114, 112, 106, 111, 111, 112, 107, 111, 120, 110, 111, 119, 112,
111, 118, 111, 119, 111, 111, 127, 111, 118, 111, 111, 119, 113, 111, 111, 119,
118, 111, 119, 113, 111, 115, 109, 114, 106, 111, 119, 112, 103, 111, 119, 114,
111, 113, 108, 111, 117, 109, 114, 110, 111, 120, 107, 111, 115, 111, 115, 111,
                        119, 112, 110, 126, 111, 118, 111, 111, 119, 116, 111,
112, 103, 115, 111, 120, 110, 111, 119, 118, 111, 111, 111, 122, 111, 111, 111,
111, 119, 114, 111, 112, 108, 111, 111, 116, 111, 111, 114, 103, 111, 111, 111,
115, 111, 111, 116, 111, 119, 118, 111, 113, 108, 112, 103, 111, 119, 111, 114,
111, 112, 106, 103, 115, 110, 111, 118, 111, 119, 111, 111, 119, 111, 118, 111,
111}; length: 301

*/
