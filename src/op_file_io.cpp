/* #include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include "basis/mcutils/parsing.h"
#include "op_file_io.h"
#include "chiral.h"

namespace chiral
{
    void read_input_to_params(MatrixParameters& params,
			      const int J0,
			      const int T0)
    {
	std::string line;
	int line_count = 0;

	// Line 1: Operator info
	{
	    ++line_count;
	    std::getline(std::cin, line);
	    
	    std::istringstream line_stream(line);
	    std::string name, order;
	    line_stream >> name >> order;
	    
	    mcutils::ParsingCheck(line_stream, line_count, line);
	    
	    std::map<std::string, ChiralOperator::Name> name_map =
		{
		    {"rc_sq", ChiralOperator::Name::charge_radius},
		    {"identity", ChiralOperator::Name::identity},
		};

	    std::map<std::string, ChiralOperator::Order> order_map =
		{
		    {"lo", ChiralOperator::Order::lo},
		    {"nlo", ChiralOperator::Order::nlo},
		    {"n2lo", ChiralOperator::Order::n2lo},
		    {"n3lo", ChiralOperator::Order::n3lo},
		    {"n4lo", ChiralOperator::Order::n4lo},
		    {"full", ChiralOperator::Order::full},
		};

	    auto search_name = name_map.find(name);
	    auto search_order = order_map.find(order);

	    if (search_name != name_map.end()
		&& search_order != order_map.end())
	    {
		params.op.name = name_map[name];
		params.op.order = order_map[order];
	    }
	    else
	    {
		mcutils::ParsingError(line_count, line,
				      "Given operator version does not exist."); 
	    }
	}

	// Line 2: Tensor and isostensor limits
	{
	    ++line_count;
	    std::getline(std::cin, line);
	    std::istringstream line_stream(line);
	    line_stream >> params.J0_min
			>> params.J0_max
			>> params.T0_min
			>> params.T0_max;
	    
	    mcutils::ParsingCheck(line_stream, line_count, line);

	    if ((0 > params.J0_min) || (params.J0_min > params.J0_max)
		|| (params.J0_max > J0))
		mcutils::ParsingError(line_count, line,
				      "Invalid J0 range specified.");

	    if ((0 > params.T0_min) || (params.T0_min > params.T0_max)
		|| (params.T0_max > T0))
		mcutils::ParsingError(line_count, line,
				      "Invalid T0 range specified.");
	}

	// Line 3: basis parameters
	{
	    ++line_count;
	    std::getline(std::cin, line);
	    std::istringstream line_stream(line);
	    line_stream >> params.Nmax
			>> params.Jmax
			>> params.hw;

	    mcutils::ParsingCheck(line_stream, line_count, line);

	    if ((params.Nmax < 0) || (params.Jmax < 0))
		mcutils::ParsingError(line_count, line,
				      "Invalid cutoff parameters");

	    if (params.hw < 0)
		mcutils::ParsingError(line_count, line,
				      "Non-physical basis scale");
	}
    }
    } */
