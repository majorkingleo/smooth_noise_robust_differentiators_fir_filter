#pragma once

#include <array>
#include <limits>
#include <stdexcept>
#include <type_traits>

/*
 * Smooth Noise Robust Differentiators Fir Filter
 *
 * Implemented according to this paper
 * http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/smooth-low-noise-differentiators/
 *
 * The Polynom is automatically calculated at compile time.
 *
 * Filter::SNRDFir::Filter<int64_t,int64_t,27*2+1> filter;
 *
 * while( ... )
 *   filter.add( adc_value );
 *
 * current_result = filter.get_result();
 *
 */

namespace exmath::Filter::SNRDFir {

namespace internal {

	template<typename T, unsigned n>
	constexpr std::array<T,n> calc_last_line_of_catalan_triangle();

	/*
	template<unsigned N> constexpr int square() { return N * N; }
	template<unsigned N> constexpr bool isodd() { return N % 2 ==1; }
	*/
	template<unsigned N> struct is_odd { static const bool value = is_odd<N-2>::value; };

	template<unsigned N> constexpr bool is_odd_v = is_odd<N>::value;

	template<> struct is_odd<1> { static const bool value = true; };

	template<> struct is_odd<0> { static const bool value = false; };

	template<typename T, T N>
	concept odds_only =
	    std::is_integral_v<T>
	    && is_odd_v<N>;


} // namespace internal

template <typename T, typename C, unsigned N>
requires internal::odds_only<unsigned, N>
class Filter
{
protected:
	std::array<T, N> input_buffer{};   // Input buffer
	std::array<C, N> coefficients = calc_coefficients();

    C       sum = 0;
    int     index = 0;
    C       default_denominator = calc_default_denominator();

public:
	/**
	 * add data without calculating
	 */
	void add(T input)
	{
		input_buffer[index] = input;
		index = (index + 1) % N;
	}

	/**
	 * Calculates the sum and stores it.
	 */
	C calculate()
	{
		C output = 0;

		/* default way to straight to accumulate
			for (int i = 0; i < N; i++) {
				output += coefficients[i] * input_buffer[(index + i) % N];
			}
		*/

		/**
		 * from outside to inside, to get accumulate data with increasing from low to high
		 * Will only matter, if you are using realy small data and float, or double as T.
		 * The center number is unused.
		 */
		for (int i = 0, j = N-1; i < N/2; i++, --j ) {
			output += coefficients[i] * input_buffer[(index + i) % N];
			output += coefficients[j] * input_buffer[(index + j) % N];
		}

		sum = output;

		return sum;
	}

	/**
	 * adds the new input value, calculates the filter and returns the devided result
	 */
	T operator()( T input )
	{
		add( input );
		calculate();
		return get_result();
	}

	/**
	 * calculate and return the devided result
	 */
	T get_result() {
		calculate();
		return sum / default_denominator;
	}

	/**
	 * choose last sum and divide it
	 */
	T get_last_result() const {
		return sum / default_denominator;
	}

	C get_default_denominator() const {
		return default_denominator;
	}

	void set_default_denominator( C dd ) {
		default_denominator = dd;
	}

	const std::array<C, N> & get_coefficients() const {
		return coefficients;
	}

	/**
	 * Tests if the given maximum input value will overflow within the calculation
	 * Throw's an exception if calculation is not possible.
	 * When result evaluated as constexpr an overflow will lead to a compiler error.
	 *
	 * Eg: constexpr auto c = filter.check_will_overflow( 4096 );
	 */
	static constexpr C check_will_it_overflow( T max_input_value )
	{
		C output = 0;
		std::array<C, N> coefficients = calc_coefficients();

		for( unsigned i = 0; i < N; ++i ) {

			C co = coefficients[i];

			C val = co * max_input_value;

			if( val > 0 ) {
				if( ( std::numeric_limits<C>::max() - val ) < val ) {
					throw std::overflow_error("Overflow error. Maximum input value too large for this datatype.");
				}
			} else {
				if( ( std::numeric_limits<C>::min() - val ) > val ) {
					throw std::overflow_error("Overflow error. Multiplication not possible with this datatype.");
				}
			}

			output += val;
		}

		return output;
	}

private:
	static constexpr std::array<C, N> calc_coefficients()
	{
		std::array<C, N> coefficients{};
		constexpr auto catalans_triangle = internal::calc_last_line_of_catalan_triangle<T,N/2>();

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

		return coefficients;
	}

	static constexpr C calc_default_denominator()
	{
		// calculate as constexpr to get an overflow error, if calculation is not possible
		constexpr C denominator = 8 * ipow<C>(4,N/2-2);
		return denominator;
	}

public:
	template <typename X>
	static constexpr X ipow(X num, unsigned int pow)
	{
		if( pow == 0 ) {
			return 1;
		}

		X mul = ipow( num, pow-1 );

		if( ( std::numeric_limits<X>::max() / mul ) < num ) {
			throw std::overflow_error("Overflow error. Multiplication not possible with this datatype.");
		}


		return num * mul;
	}
};

namespace internal {

/*
 * Calcualates the last line of the catalan triangle
 * https://oeis.org/A039598
 * https://oeis.org/A039598/b039598.txt
 *
 * (Sage)
 *  D = [0]*(n+2); D[1] = 1
 *  b = True; h = 1
 *  for i in range(2*n) :
 *      if b :
 *        for k in range(h, 0, -1) : D[k] += D[k-1]
 *          h += 1
 *      else :
 *          for k in range(1, h, 1) : D[k] += D[k+1]
 *      b = not b
 *      if b : print([D[z] for z in (1..h-1) ])
 */

template<typename T, unsigned n>
constexpr std::array<T,n> calc_last_line_of_catalan_triangle()
{
  std::array<T,n+2> D{0};

  bool b = true;
  int h = 1;
  D.at(1) = 1;

  for( int i = 0; i < 2 * int(n); ++i ) {
	if( b ) {
	  for( int k = h; k > 0; --k ) {

		if( std::numeric_limits<T>::max() - D.at(k-1) < D.at(k) ) {
		  throw std::overflow_error("Overflow error. Coefficient not possible with this datatype.");
		}

		D.at(k) += D.at(k-1);
	  }
	  h += 1;
	} else {
	  for( int k = 1; k < h; ++k ) {

		if( std::numeric_limits<T>::max() - D.at(k+1) < D.at(k) ) {
		  throw std::overflow_error("Overflow error. Coefficient not possible with this datatype.");
		}

		D.at(k) += D.at(k+1);
	  }
	}

	b = !b;
  }

  std::array<T,n> ret{0};

  for( int z = 1,i = 0; z < h; ++z, ++i ) {
	ret.at(i) = D.at(z);
  }

  return ret;
}

} // namespace internal

} // namespace Filter::SNRD
