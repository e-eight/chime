#ifndef THREEDHO_H
#define THREEDHO_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>
#include "constants.h"

namespace ho
{
  // Oscillator parameter for relative & relative-cm cases.
  struct OscillatorParameter
  {
    OscillatorParameter()
      : oscillator_energy(0) {}
    OscillatorParameter(double _energy)
      : oscillator_energy(_energy) {}
    ~OscillatorParameter() = default;

    double relative() const
    {
      auto b = constants::hbarc;
      b /= std::sqrt(constants::reduced_nucleon_mass_MeV * oscillator_energy);
      return b;
    }
    double cm() const
    {
      auto b = constants::hbarc;
      b /= std::sqrt(2 * constants::nucleon_mass_MeV * oscillator_energy);
      return b;
    }
  private:
    double oscillator_energy; // in MeV
  };

  // Coordinate space normalization constant.
  static inline double CoordinateSpaceNorm(const std::size_t n,
                                           const std::size_t l,
                                           const double b)
  {
    auto result = (std::log(2) + gsl_sf_lngamma(n+1)
                   - 3 * std::log(b) - gsl_sf_lngamma(n+l+1.5));
    return std::sqrt(std::exp(result));
  }

  // Momentum space normalization constant.
  static inline double MomentumSpaceNorm(const std::size_t n,
                                         const std::size_t l,
                                         const double b)
  {
    auto result = (std::log(2) + gsl_sf_lngamma(n+1)
                   + 3 * std::log(b) - gsl_sf_lngamma(n+l+1.5));
    return std::sqrt(std::exp(result));
  }


  // 3DHO wave function using recursion.
  //
  // In coordinate space:
  // Ψ(n, l, r, b) = (A(n-1, l+1/2) (r/b)^2 + B(n-1, l+1/2)) Ψ(n-1, l, r, b)
  //                   - C(n-1, l+1/2) Ψ(n-2, l, r, b),
  // Ψ(0, l, r, b) = \sqrt{2/b^3} exp(-r^2/2b^2) (r/b)^l \sqrt{1/Γ(3/2+l)},
  // Ψ(1, l, r, b) = \sqrt{2/b^3} exp(-r^2/2b^2) (r/b)^l (3/2 - r^2/b^2 + l) \sqrt{1/Γ(5/2+l)}.
  //
  // In momentum space:
  // Ψ(n, l, p, b) = (A(n-1, l+1/2) (bp)^2 + B(n-1, l+1/2)) Ψ(n-1, l, p, b)
  //                   - C(n-1, l+1/2) Ψ(n-2, l, p, b),
  // Ψ(0, l, r, b) = \sqrt{2b^3} exp(-b^2p^2/2) (bp)^l \sqrt{1/Γ(3/2+l)},
  // Ψ(1, l, r, b) = \sqrt{2b^3} exp(-b^2p^2/2) (bp)^l (3/2 - b^2p^2 + l) \sqrt{1/Γ(5/2+l)}.
  //
  // The coefficients are:
  // A(n, a) = - 1 / \sqrt{(1+n)(1+a+n)}
  // B(n, a) = (1+a+2n) / \sqrt{(1+n)(1+a+n)}
  // C(n, a) = n(a+n) / \sqrt{n(1+n)(a+n)(1+a+n)}

  // double AW(int n, double a)
//   {
//     return -1.0 / std::sqrt((n + 1) * (n + a + 1));
//   }
//   double BW(int n, double a)
//   {
//     return (2*n + a + 1) / std::sqrt((n + 1) * (n + a + 1));
//   }
//   double CW(int n, double a)
//   {
//     return std::sqrt(n * (n + a) / ((n + 1) * (n + a + 1)));
//   }

//   // Calculate the value of the 3DHO basis function at a given point.
//   double WaveFunctionValue(int n, int l, double b, double x)
//   {
//     if (n == 0)
//       {
//         double value = std::sqrt(2.0 / (cube(b) * gsl_sf_gamma(l + 1.5)));
//         value *= std::exp(-0.5 * square(x / b)) * std::pow(x / b, l);
//         return value;
//       }
//     if (n == 1)
//       {
//         double value = std::sqrt(2.0 / (cube(b) * gsl_sf_gamma(l + 2.5)));
//         value *= (std::exp(-0.5 * square(x / b)) * std::pow(x / b, l)
//                   * (l + 1.5 - square(x / b)));
//         return value;
//       }

//     std::vector<double> table(n + 1);
//     table[0] = WaveFunctionValue(0, l, b, x);
//     table[1] = WaveFunctionValue(1, l, b, x);
//     for (int i = 2; i <= n; ++i)
//       {
//         table[i] = ((AW(i-1, l+0.5) * square(x / b) + BW(i-1, l+0.5)) * table[i-1]
//                     - CW(i-1, l+0.5) * table[i-2]);
//       }
//     return table.back();
//   }

//   // Calculates the value of the 3DHO basis function at the given vector of points.
//   std::vector<double> WaveFunctionValue(int n, int l, double b, std::vector<double> xs)
//     {
//       int m = xs.size();
//       std::vector<double> values(m);
//       std::transform(xs.begin(), xs.end(), values.begin(),
//                      [=](double x)->double{return WaveFunctionValue(n, l, b, x);});
//       return values;
//     }
// }
}

#endif
