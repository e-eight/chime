#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include "basis/mcutils/parsing.h"
#include "op_file_io.h"
#include "chiral.h"

namespace chiral
{

    MatrixParameters input_to_params()
    {
	MatrixParameters params;
	std::string line;
	int line_count = 0;
  
	{
	    ++line_count;
	    std::getline(std::cin, line);
	    std::istringstream line_stream(line);
	
	    line_stream >> params.name >> params.order;
	  
	    params.op = ChiralOperator::create_operator(params.name);
	    std::map<std::string, Order> orders =
		{
		    {"lo", Order::lo},
		    {"nlo", Order::nlo},
		    {"n2lo", Order::n2lo},
		    {"n3lo", Order::n3lo},
		    {"n4lo", Order::n4lo},
		    {"full", Order::full},
		};
	    params.op->order = orders[params.order];
	}
	
	{
	    ++line_count;
	    std::getline(std::cin, line);
	    std::istringstream line_stream(line);

	    line_stream >> params.J0_min
			>> params.J0_max
			>> params.T0_min
			>> params.T0_max;

	}

	{
	    ++line_count;
	    std::getline(std::cin, line);
	    std::istringstream line_stream(line);

	    line_stream >> params.Nmax
			>> params.Jmax
			>> params.hw;
	}
	
	return params;
    }

}
