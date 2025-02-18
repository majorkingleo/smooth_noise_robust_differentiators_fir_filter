#include "catalans_triangle.h"
#include <iostream>
#include <ColoredOutput.h>
#include <OutDebug.h>
#include <arg.h>
#include <file_option.h>
#include <stderr_exception.h>
#include <fstream>
#include "SNRDFir.hpp"

using namespace Tools;
using namespace exmath;

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

			Filter::SNRDFir::Filter<int64_t,int64_t,27*2+1> filter;
			filter.set_default_denominator(filter.get_default_denominator()/ 256);

			constexpr auto c = filter.check_will_it_overflow( 4096 );
			/*
			if( filter.check_will_overflow(1) ) {
				throw std::out_of_range( "value to large" );
			}*/

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

			Filter::SNRDFir::Filter<float,float,63*2+1> filter;
			filter.set_default_denominator(filter.get_default_denominator()/ 256.0);
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

			Filter::SNRDFir::Filter<double,double,397*2+1> filter;
			filter.set_default_denominator(filter.get_default_denominator()/ 256.0);
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

