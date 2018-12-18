#include <iostream>
#include "operator_rc_sq.h"
// #include "utility.h"

int main()
{
    for (int ml = -2; ml <= 2; ++ml)
    {
	std::cout << gsl_sf_coupling_3j(4, 2, 2, 2*ml, -2*ml, 0) << std::endl;
    }
}

