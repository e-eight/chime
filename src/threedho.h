#ifndef THREEDHO_H
#define THREEDHO_H

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

namespace threedho
{
  inline double CoordinateSpaceNorm(int n, int l, double b)
  {
    double numerator = M_LN2 + gsl_sf_lngamma(n + 1);
    double denominator = 3 * std::log(b) + gsl_sf_lngamma(n + l + 1.5);
    return std::exp(0.5 * (numerator - denominator));
  }

  inline double MomentumSpaceNorm(int n, int l, double b)
  {
    double numerator = M_LN2 + gsl_sf_lngamma(n + 1) + 3 * std::log(b);
    double denominator = gsl_sf_lngamma(n + l + 1.5);
    return std::exp(0.5 * (numerator - denominator));
  }
}

#endif
