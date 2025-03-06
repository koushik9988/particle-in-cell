#ifndef _Fourier_transform_
#define _Fourier_transform_

#include <vector>
#include <cmath>
#include "slap.h"

namespace dft
{
    constexpr double PI = 3.14159265358979323846;

    struct dft_result
    {
        vec<double> freq;
        std::vector<vec<double>> dft_amplitude;
        vec<double> magnitude;            
    };

    dft_result dft_1d(vec<double> &f, double dx);
}

#endif