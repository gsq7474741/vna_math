//   fft.cpp - impelementation of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

//   http://www.librow.com/articles/article-10

//   Reworked from original from LIBROW (info above) to Arduino Lib in c++

//   Include declaration file
#include <fast4ier.h>


//   FORWARD FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
bool Fast4::FFT(const complex *const Input, complex *const Output, const unsigned int n)
{
	//   Check input parameters
	if (!Input || !Output || n < 1 || n & (n - 1))
		return false;
	//   Initialize data
	Rearrange(Input, Output, n);
	//   Call FFT implementation
	Perform(Output, n);
	//   Succeeded
	return true;
}

//   FORWARD FOURIER TRANSFORM, INPLACE VERSION
//     Data - both input data and output
//     N    - length of input data
bool Fast4::FFT(complex *const Data, const unsigned int n)
{
	//   Check input parameters
	if (!Data || n < 1 || n & (n - 1))
		return false;
	//   Rearrange
	Rearrange(Data, n);
	//   Call FFT implementation
	Perform(Data, n);
	//   Succeeded
	return true;
}

//   INVERSE FOURIER TRANSFORM
//     Input  - input data
//     Output - transform result
//     N      - length of both input data and result
//     Scale  - if to scale result
bool Fast4::IFFT(const complex *const Input, complex *const Output, const unsigned int n, const bool Scale /* = true */)
{
	//   Check input parameters
	if (!Input || !Output || n < 1 || n & (n - 1))
		return false;
	//   Initialize data
	Rearrange(Input, Output, n);
	//   Call FFT implementation
	Perform(Output, n, true);
	//   Scale if necessary
	if (Scale)
		Fast4::Scale(Output, n);
	//   Succeeded
	return true;
}

//   INVERSE FOURIER TRANSFORM, INPLACE VERSION
//     Data  - both input data and output
//     N     - length of both input data and result
//     Scale - if to scale result
bool Fast4::IFFT(complex *const Data, const unsigned int n, const bool Scale /* = true */)
{
	//   Check input parameters
	if (!Data || n < 1 || n & (n - 1))
		return false;
	//   Rearrange
	Rearrange(Data, n);
	//   Call FFT implementation
	Perform(Data, n, true);
	//   Scale if necessary
	if (Scale)
		Fast4::Scale(Data, n);
	//   Succeeded
	return true;
}

//   Rearrange function
void Fast4::Rearrange(const complex *const Input, complex *const Output, const unsigned int n)
{
	//   Data entry position
	unsigned int Target = 0;
	//   Process all positions of input signal
	for (unsigned int Position = 0; Position < n; ++Position)
	{
		//  Set data entry
		Output[Target] = Input[Position];
		//   Bit mask
		unsigned int Mask = n;
		//   While bit is set
		while (Target & (Mask >>= 1))
			//   Drop bit
			Target &= ~Mask;
		//   The current bit is 0 - set it
		Target |= Mask;
	}
}

//   Inplace version of rearrange function
void Fast4::Rearrange(complex *const Data, const unsigned int n)
{
	//   Swap position
	unsigned int Target = 0;
	//   Process all positions of input signal
	for (unsigned int Position = 0; Position < n; ++Position)
	{
		//   Only for not yet swapped entries
		if (Target > Position)
		{
			//   Swap entries
			const complex Temp(Data[Target]);
			Data[Target] = Data[Position];
			Data[Position] = Temp;
		}
		//   Bit mask
		unsigned int Mask = n;
		//   While bit is set
		while (Target & (Mask >>= 1))
			//   Drop bit
			Target &= ~Mask;
		//   The current bit is 0 - set it
		Target |= Mask;
	}
}

//   FFT implementation
void Fast4::Perform(complex *const Data, const unsigned int n, const bool Inverse /* = false */)
{
	const FLT pi = Inverse ? 3.14159265358979323846 : -3.14159265358979323846;
	//   Iteration through dyads, quadruples, octads and so on...
	for (unsigned int Step = 1; Step < n; Step <<= 1)
	{
		//   Jump to the next entry of the same transform factor
		const unsigned int Jump = Step << 1;
		//   Angle increment
		const FLT delta = pi / FLT(Step);
		//   Auxiliary sin(delta / 2)
		const FLT Sine = sin(delta * .5);
		//   Multiplier for trigonometric recurrence
		const complex Multiplier(-2. * Sine * Sine, sin(delta));
		//   Start value for transform factor, fi = 0
		complex Factor(1.);
		//   Iteration through groups of different transform factor
		for (unsigned int Group = 0; Group < Step; ++Group)
		{
			//   Iteration within group 
			for (unsigned int Pair = Group; Pair < n; Pair += Jump)
			{
				//   Match position
				const unsigned int Match = Pair + Step;
				//   Second term of two-point transform
				const complex Product(Factor * Data[Match]);
				//   Transform for fi + pi
				Data[Match] = Data[Pair] - Product;
				//   Transform for fi
				Data[Pair] += Product;
			}
			//   Successive transform factor via trigonometric recurrence
			Factor = Multiplier * Factor + Factor;
		}
	}
}

//   Scaling of inverse FFT result
void Fast4::Scale(complex *const Data, const unsigned int n)
{
	const FLT Factor = 1. / FLT(n);
	//   Scale all data entries
	for (unsigned int Position = 0; Position < n; ++Position)
		Data[Position] *= Factor;
}
