//   fft.h - declaration of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

//   http://www.librow.com/articles/article-10

//   Reworked from original from LIBROW (info above) to Arduino Lib in c++

#ifndef fast4ier_h
#define fast4ier_h
//#include "Arduino.h"
//   Include complex numbers header
#include <complex.h>
//   Include math library
#include "math.h"
class Fast4
{
public:
	//   FORWARD FOURIER TRANSFORM
	//     Input  - input data
	//     Output - transform result
	//     N      - length of both input data and result
	static bool FFT(const complex *const Input, complex *const Output, const unsigned int n);

	//   FORWARD FOURIER TRANSFORM, INPLACE VERSION
	//     Data - both input data and output
	//     N    - length of input data
	static bool FFT(complex *const Data, const unsigned int n);

	//   INVERSE FOURIER TRANSFORM
	//     Input  - input data
	//     Output - transform result
	//     N      - length of both input data and result
	//     Scale  - if to scale result
	static bool IFFT(const complex *const Input, complex *const Output, const unsigned int n, const bool Scale = true);

	//   INVERSE FOURIER TRANSFORM, INPLACE VERSION
	//     Data  - both input data and output
	//     N     - length of both input data and result
	//     Scale - if to scale result
	static bool IFFT(complex *const Data, const unsigned int n, const bool Scale = true);

private:
	//   Rearrange function and its inplace version
	static void Rearrange(const complex *const Input, complex *const Output, const unsigned int n);
	static void Rearrange(complex *const Data, const unsigned int n);

	//   FFT implementation
	static void Perform(complex *const Data, const unsigned int n, const bool Inverse = false);

	//   Scaling of inverse FFT result
	static void Scale(complex *const Data, const unsigned int n);
};

#endif
