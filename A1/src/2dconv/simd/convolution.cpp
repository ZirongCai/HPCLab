#include <iostream>
#include <chrono>
#include <immintrin.h>

#include "compute_kernels.h"

using namespace std::chrono;

namespace conv
{
    void apply_simd_intrinsics(ImageRGB &input,
                               Filter &filter,
                               ImageRGB &output)
    {
        if (input != output)
        {
            std::cout << "ERROR: input and output dimensions do not match" << std::endl;
            REPORT_ERR;
            exit(EXIT_FAILURE);
        }

        float *i_channels[ImageRGB::num_channels] __attribute__((aligned(64))) = {input.get_red(), input.get_green(), input.get_blue()};
        float *o_channels[ImageRGB::num_channels] = {output.get_red(), output.get_green(), output.get_blue()};
        const float *filter_data = filter.get_data();

        const int k_half_width = filter.get_half_width();
        const int k_half_height = filter.get_half_height();
        const int k_total_size = filter.get_total_size();

        const int input_size = input.get_width() * input.get_height();

        const int register_length = 8;

        CommandLineSettings *settings = CommandLineSettings::get_settings();

        auto start = steady_clock::now();

        //
        // TODO: implement your code
        //
        int image_end = input.get_width() - input.get_width() % register_length;

        for (size_t repeat = 0; repeat < settings->get_num_repeats(); ++repeat)
        {
            for (size_t channel = 0; channel < ImageRGB::num_channels; ++channel)
            {
                // for (size_t index = 0; index < image_end; index += 8)
                // {
                for (size_t y = 0; y < input.get_height(); ++y)
                {
                    for (size_t x = 0; x < image_end; ++x)
                    {
                        // int y = index / input.get_width();
                        // int x = index % input.get_width();
                        //caculate 8 output at once
                        __m256 C; // filter coefficinet vector
                        __m256 I; // image data vector
                        __m256 result = _mm256_setzero_ps();
                        for (int yy = -k_half_height; yy <= k_half_height; ++yy)
                        {
                            for (int xx = -k_half_width; xx <= k_half_width; ++xx)
                            {
                                C = _mm256_set1_ps(filter_data[filter.get_rel_index(xx, yy)]);
                                I = _mm256_loadu_ps(&i_channels[channel][input.get_index(x + xx, y + yy)]);
                                result = _mm256_fmadd_ps(C, I, result);
                            }
                        }

                        _mm256_storeu_ps(&o_channels[channel][output.get_index(x, y)], result);
                        // }
                    }
                    //handle Remainder
                    for (size_t remainder = image_end; remainder < input.get_width(); ++remainder)
                    {
                        // int y = remainder / input.get_width();
                        // int x = remainder % input.get_width();
                        float value = 0.0f;
                        for (int yy = -k_half_height; yy <= k_half_height; ++yy)
                        {
                            for (int xx = -k_half_width; xx <= k_half_width; ++xx)
                            {

                                value += (i_channels[channel][input.get_index(remainder + xx, y + yy)] * filter_data[filter.get_rel_index(xx, yy)]);
                            }
                        }
                        o_channels[channel][output.get_index(remainder, y)] = value;
                    }
                }
            }
        }

        auto end = steady_clock::now();
        duration<double> elapse_time = end - start;
        std::cout << "time spent: " << elapse_time.count() << ", sec" << std::endl;

        // float test_value = 0;
        // for (size_t y = 0; y < 10; ++y)
        // {
        //     for (size_t x = 0; x < 10; ++x)
        //     {
        //         test_value += o_channels[0][output.get_index(x, y)];
        //     }
        // }
        // std::cout << "Default Output sum = " << test_value << std::endl;
    }

    void apply_asm(ImageRGB &input, Filter &filter, ImageRGB &output)
    {
        std::cout << "convolution with asm has not been implemented" << std::endl;
    }
}
