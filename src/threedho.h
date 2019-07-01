#ifndef THREEDHO_H
#define THREEDHO_H

#include <cmath>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

namespace threedho
{
  inline double CoordinateSpaceNorm(const int& n, const int& l, const double& b)
  {
    double numerator = M_LN2 + gsl_sf_lngamma(n + 1);
    double denominator = 3 * std::log(b) + gsl_sf_lngamma(n + l + 1.5);
    return std::exp(0.5 * (numerator - denominator));
  }

  inline double MomentumSpaceNorm(const int& n, const int& l, const double& b)
  {
    double numerator = M_LN2 + gsl_sf_lngamma(n + 1) + 3 * std::log(b);
    double denominator = gsl_sf_lngamma(n + l + 1.5);
    return std::exp(0.5 * (numerator - denominator));
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

  inline double CoefficientA(const int& n, const double& a)
  {
    return (-1 / std::sqrt((1 + n) * (1 + a + n)));
  }

  inline double CoefficientB(const int& n, const double& a)
  {
    return ((1 + a + 2*n) / std::sqrt((1 + n) * (1 + a + n)));
  }

  inline double CoefficientC(const int& n, const double& a)
  {
    return (n * (a + n) / std::sqrt(n * (1 + n) * (a + n) * (1 + a + n)));
  }

  enum struct Space { configuration, momentum };
  struct WaveFunction
  {
    WaveFunction(): n(0), l(0), b(1.0), space(Space::configuration) {}
    WaveFunction(int n, int l, double b, Space space)
      : n(n), l(l), b(b), space(space) {}
    ~WaveFunction(){}

    double operator()(double x)
    {
      auto bx = (space == Space::configuration ? x / b : b * x);
      auto prefactor = (space == Space::configuration
                        ? std::sqrt(2 / (b * b * b))
                        : std::sqrt(2 * b * b * b));
      std::vector<double> table(n+1);
      table[0] = (prefactor * std::exp(bx * bx / 2) * std::pow(bx, l)
                  * std::sqrt(1 / std::tgamma(l + 1.5)));
      table[1] = (prefactor * std::exp(bx * bx / 2) * std::pow(bx, l)
                  * (1.5 - bx + l) * std::sqrt(1 / std::tgamma(l + 2.5)));
      for (int i = 2; i <= n; ++i)
        {
          table[i] = ((CoefficientA(i - 2, l + 0.5) * bx * bx
                      + CoefficientB(i - 1, l + 0.5)) * table[i - 1]
                      - CoefficientC(i - 1, l + 0.5) * table[i - 2]);
        }
      return table[n];
    }

  private:
    int n;
    int l;
    double b;
    Space space;
  };
}

#endif
