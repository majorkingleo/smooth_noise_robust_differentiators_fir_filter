#include "catalans_triangle.h"
#include <iostream>

template<class T>
static void dump_array( const T & co )
{
	std::cout << "[";
	for( unsigned i = 0; i < co.size(); ++i  ) {
		if( i > 0 ) {
			std::cout << ", ";
		}

		std::cout << co[i];
	}

	std::cout << "]\n";
}

template <typename T, typename C, unsigned N>
class Filter
{
protected:
	std::array<T, N> input_buffer{};   // Input buffer
	std::array<C, N> coefficients{};

    C sum = 0;
    int index = 0;

public:
	constexpr Filter()
	{
		init_coefficients();
	}

	constexpr T getOutput() const { return sum; }

	constexpr T operator()(T input)
	{
		// Convert input to float for processing
		input_buffer[index] = input;

		C output = 0;
		for (int i = 0; i < N; i++) {
			sum += coefficients[i] * input_buffer[(index + i) % N];
		}

		sum = output;
		index = (index + 1) % N;
		return sum;
	}

	const std::array<C, N> & get_coefficients() const {
		return coefficients;
	}

private:
	constexpr void init_coefficients()
	{
		constexpr auto catalans_triangle = calc_last_line_of_catalan_triangle<T,N/2>();

		// fill reverse and negative
		unsigned i = 0;
		for( auto it = catalans_triangle.rbegin(); it < catalans_triangle.rend(); ++it, ++i ) {
			coefficients[i] = *it * -1;
		}

		// the middle one is unused
		coefficients[i] = 0;
		++i;

		for( auto it = catalans_triangle.begin(); it < catalans_triangle.end(); ++it, ++i ) {
			coefficients[i] = *it;
		}
	}
};

template<class T>
static void dump_coefficients( const T & filter )
{
	dump_array( filter.get_coefficients() );
}

int main( int arg, char **argv )
{
	try {
		Filter<int64_t,int64_t,27*2+1> f;

		dump_coefficients(f);

		std::cout << "(uint64_t(1) << 29) = " << (uint64_t(1) << 29) << std::endl;

	} catch( const std::exception & error ) {
		std::cerr << "Error: " << error.what() << std::endl;
	}

	return 0;
}

