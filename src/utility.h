#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <cmath>
#include <gsl/gsl_sf_coupling.h>

namespace util
{    
    inline int delta(int a, int b) { return a == b; };

    inline double clebsch(int ja, int jb, int j, int mja, int mjb, int mj)
    {
	return (pow(-1, ja - jb + mj) * sqrt(2 * j + 1)
		* gsl_sf_coupling_3j(2 * ja, 2 * jb, 2 * j,
				     2 * mja, 2 * mjb, -2 * mj));
    }
}

#endif
