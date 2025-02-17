#include "catalans_triangle.h"
#include <iostream>
#include <ColoredOutput.h>
#include <OutDebug.h>
#include <arg.h>
#include <file_option.h>
#include <stderr_exception.h>
#include <fstream>

using namespace Tools;

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

    C       sum = 0;
    int     index = 0;
    C       default_devide = 0;
    bool    start_shifting = false;

public:
	constexpr Filter()
	{
		init_coefficients();
	}

	constexpr T get_result() const {
		// if( sum < 0 ) {
		//	CPPDEBUG( format( "sum: %d / %d", sum, default_devide ));
		//}

		return sum / default_devide;
	}

	constexpr T operator()(T input)
	{
		// Convert input to float for processing
		input_buffer[index] = input;

		C output = 0;
		/*
		for (int i = 0; i < N; i++) {
			output += coefficients[i] * input_buffer[(index + i) % N];
		}*/

		for (int i = 0, j = N-1; i < N/2; i++, --j ) {
			//CPPDEBUG( format( "i: %d j: %d", i, j ) );
			output += coefficients[i] * input_buffer[(index + i) % N];
			output += coefficients[j] * input_buffer[(index + j) % N];
		}


		sum = output;
		index = (index + 1) % N;

		return sum;
	}

	const std::array<C, N> & get_coefficients() const {
		return coefficients;
	}

	C get_default_devide() const {
		return default_devide;
	}

	void set_default_devide( C dd ) {
		default_devide = dd;
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

		// calculate as constexpr to get an overflow error, if calculation is not possible
		constexpr C devider = 8 * ipow<C>(4,N/2-2);
		default_devide = devider;
	}

public:
	template <typename X>
	static constexpr X ipow(X num, unsigned int pow)
	{
		if( pow == 0 ) {
			return 1;
		}

		X mul = ipow( num, pow-1 );

		if( std::numeric_limits<X>::max() / mul < num ) {
			throw std::overflow_error("Overflow error. Multiplication not possible with this datatype.");
		}


		return num * mul;
	}
};

template<class T>
static void dump_coefficients( const T & filter )
{
	dump_array( filter.get_coefficients() );
}

static int32_t get_as_12bit_adc( float volt )
{
	const unsigned ADC_MAX12 = 0x1000;
	int32_t adc_value = volt * double(ADC_MAX12) / 3.3;

	return adc_value;
}

static float get_12bit_adc_as_volt( int32_t adc_value )
{
	const unsigned ADC_MAX12 = 0x1000;
	float volt = adc_value * 3.3 / double(ADC_MAX12);

	return volt;
}


int main( int argc, char **argv )
{
	try {
		ColoredOutput co;

		Arg::Arg arg( argc, argv );
		arg.addPrefix( "-" );
		arg.addPrefix( "--" );

		Arg::OptionChain oc_info;
		arg.addChainR(&oc_info);
		oc_info.setMinMatch(1);
		oc_info.setContinueOnMatch( false );
		oc_info.setContinueOnFail( true );

		Arg::FlagOption o_help( "help" );
		o_help.setDescription( "Show this page" );
		oc_info.addOptionR( &o_help );

		Arg::FlagOption o_debug("d");
		o_debug.addName( "debug" );
		o_debug.setDescription("print debugging messages");
		o_debug.setRequired(false);
		arg.addOptionR( &o_debug );

		Arg::FlagOption o_fir1("fir1");
		o_fir1.setDescription("Ulfs suggested FIR filter");
		o_fir1.setRequired(false);
		arg.addOptionR( &o_fir1 );

		Arg::FlagOption o_fir2("fir2");
		o_fir2.setDescription("FIR filter with float 67 cofficients");
		o_fir2.setRequired(false);
		arg.addOptionR( &o_fir2 );

		Arg::FlagOption o_fir3("fir3");
		o_fir3.setDescription("FIR filter with double 67 cofficients");
		o_fir3.setRequired(false);
		arg.addOptionR( &o_fir3 );

		Arg::EmptyFileOption o_file;
		o_file.setDescription("input file");
		o_file.setRequired(true);
		arg.addOptionR( &o_file );

		if( !arg.parse() )
		{
			std::cout << arg.getHelp(5,20,30, 80 ) << std::endl;
			return 1;
		}

		if( o_debug.getState() )
		{
			Tools::x_debug = new OutDebug();
		}

		if( o_help.getState() ) {
			std::cout << arg.getHelp(5,20,30, 80 ) << std::endl;
			return 1;
		}

		if( o_fir1.getState() ) {

			std::ifstream in( o_file.getValues()->at(0) );

			if( !in ) {
				throw STDERR_EXCEPTION( Tools::format( "cannot open file %s", o_file.getValues()->at(0) ) );
			}

			Filter<int64_t,int64_t,27*2+1> filter;
			filter.set_default_devide(filter.get_default_devide()/ 256);
			//dump_coefficients(filter);

			while( !in.eof() ) {

				float f_in = 0;
				in >> f_in;

				uint32_t adc = get_as_12bit_adc( f_in );
				//std::cout << adc << std::endl;

				// std::cout << filter(adc) << std::endl;

				filter(adc);
				int64_t filtered_adc = filter.get_result();

				std::cout << get_12bit_adc_as_volt(filtered_adc) << std::endl;
				//std::cout << filtered_adc << std::endl;
			}
		}
		else if( o_fir2.getState() ) {

			std::ifstream in( o_file.getValues()->at(0) );

			if( !in ) {
				throw STDERR_EXCEPTION( Tools::format( "cannot open file %s", o_file.getValues()->at(0) ) );
			}

			Filter<float,float,63*2+1> filter;
			filter.set_default_devide(filter.get_default_devide()/ 256.0);
			//dump_coefficients(filter);

			while( !in.eof() ) {

				float f_in = 0;
				in >> f_in;

				filter(f_in);

				std::cout << filter.get_result() << std::endl;
			}

		}
		else if( o_fir3.getState() ) {

			std::ifstream in( o_file.getValues()->at(0) );

			if( !in ) {
				throw STDERR_EXCEPTION( Tools::format( "cannot open file %s", o_file.getValues()->at(0) ) );
			}

			Filter<double,double,397*2+1> filter;
			filter.set_default_devide(filter.get_default_devide()/ 256.0);
			//dump_coefficients(filter);

			while( !in.eof() ) {

				float f_in = 0;
				in >> f_in;

				filter(f_in);

				std::cout << filter.get_result() << std::endl;
			}

		}

	} catch( const std::exception & error ) {
		std::cerr << "Error: " << error.what() << std::endl;
	} catch( ... ) {
		std::cerr << "Unknown Error\n";
	}


	return 0;
}

