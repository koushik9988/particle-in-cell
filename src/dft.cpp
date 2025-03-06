#include "dft.h"

namespace dft
{
    // 1D Discrete Fourier Transform
    dft_result dft_1d(vec<double> &f, double dx)
    {
        int N = f.getSize();
        dft_result result;
        result.freq = vec<double>(N);
        result.dft_amplitude.resize(N);
        result.magnitude = vec<double>(N);

        for (int i = 0; i < N; i++)
        {
            result.dft_amplitude[i] = vec<double>(2); // Allocate space for real and imaginary parts
        }

        // Compute FT
        vec<double> ft_kernel(2);
        vec<double> sum(2);
        for (int k = 0; k < N; k++)
        {
            sum = 0; // Reset sum for each frequency component
            for (int n = 0; n < N; n++)
            {
                double theta = -2.0 * PI * k * n / static_cast<double>(N);
                ft_kernel(0) = cos(theta);
                ft_kernel(1) = sin(theta);
                sum += ft_kernel * f(n); // Apply Fourier kernel 
            }
            result.dft_amplitude[k] = sum; 
            result.magnitude(k) = sum.norm(); 
            result.freq(k) = k / (N * dx); 
        }
        return result;
    }
}