#ifndef SPECFUNC_H
#define SPECFUNC_H

#include <gsl/gsl_sf.h>

namespace sf
{
  double Laguerre(const int n, const double a, const double x)
  {
    switch (n)
      {
      case 1:
        return gsl_sf_laguerre_1(a, x);
      case 2:
        return gsl_sf_laguerre_2(a, x);
      case 3:
        return gsl_sf_laguerre_3(a, x);
      default:
        return gsl_sf_laguerre_n(n, a, x);
    }
  }
}

#endif
