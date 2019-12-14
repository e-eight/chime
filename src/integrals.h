#ifndef INTEGRALS_H
#define INTEGRALS_H

/*
 * Evaluates some of the common integrals that appear in the matrix element
 * evaluation of most operators. The radial integrals are evaluated using
 * adaptive Gauss-Kronrod quadrature from GSL.
 */

#include "quadrature.h"
#include "constants.h"
#include "specfunc.h"
#include "utility.h"
#include "threedho.h"

namespace quadrature
{
  // Lenpic coordinate space regulator.
  static inline double LenpicSemiLocalRegulator(const double r, const double R)
  {
    return std::pow(1 - std::exp(-square(r / R)), 6);
  }
  static inline double ApplyRegulator(const bool regularize,
                                      const double r,
                                      const double R)
  {
    if (regularize)
      return LenpicSemiLocalRegulator(r, R);
    else
      return 1;
  }

  // Yukawa and associated functions that appear in many operator matrix elements.
  static inline double YPi(const double r, const double mpi)
  {
    double rpi = mpi * r;
    return std::exp(-rpi) / rpi;
  }
  static inline double ZPi(const double r, const double mpi)
  {
    return 1 + mpi * r;
  }
  static inline double TPi(const double r, const double mpi)
  {
    return 2 * mpi * r - 1;
  }
  static inline double WPi(const double r, const double mpi)
  {
    double rpi = mpi * r;
    return 1 + (3 / rpi) + (3 / square(rpi));
  }

  // Integrates m_π r, sandwiched between coordinate space HO radial
  // basis states, without the regulator. This integral is always associated with
  // ⟨l'||C_1(\hat{r})||l⟩, which only allows l' = l ± 1.
  static double IntegralMPiR(const int& np,
                             const int& lp,
                             const int& n,
                             const int& l,
                             const double& b,
                             const double& mpi)
  {
    if (lp == l + 1)
      return (mpi * b) * (std::sqrt(n + l + 1.5) * (np == n) - std::sqrt(n) * (np + 1 == n));
    if (lp + 1 == l)
      return IntegralMPiR(n, l, np, lp, b, mpi);
    return 0;
  }

  static double IntegralYPiR(const int& np,
                             const int& lp,
                             const int& n,
                             const int& l,
                             const double& b,
                             const double& mpi,
                             const bool& regularize,
                             const double& R)
  {
    auto integrand =
      [&](double r)
      {
        return (ApplyRegulator(regularize, r, R)
                * YPi(r, mpi)
                * std::exp(-square(r / b))
                * std::pow(r / b, l + lp + 2)
                * sf::Laguerre(n, l + 0.5, square(r / b))
                * sf::Laguerre(np, lp + 0.5, square(r / b)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm_product = (ho::CoordinateSpaceNorm(n, l, 1)
                           * ho::CoordinateSpaceNorm(np, lp, 1) / b);
    integral *= norm_product;
    return integral;
  }

  // Integrates m_π r Y_π(r), sandwiched between coordinate space HO radial
  // basis states.
  static double IntegralMPiRYPiR(const int& np,
                                 const int& lp,
                                 const int& n,
                                 const int& l,
                                 const double& b,
                                 const double& mpi,
                                 const bool& regularize,
                                 const double& R)
  {
    auto integrand =
      [&](double r)
      {
        return (ApplyRegulator(regularize, r, R)
                * std::exp(-square(r / b) - (mpi * r))
                * std::pow(r / b, l + lp + 2)
                * sf::Laguerre(n, l + 0.5, square(r / b))
                * sf::Laguerre(np, lp + 0.5, square(r / b)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm_product = (ho::CoordinateSpaceNorm(n, l, 1)
                           * ho::CoordinateSpaceNorm(np, lp, 1) / b);
    integral *= norm_product;
    return integral;
  }

  // Integrates W_π(r) Y_π(r), sandwiched between coordinate space HO
  // radial basis states.
  static double IntegralWPiYPiR(const int& np,
                                const int& lp,
                                const int& n,
                                const int& l,
                                const double& b,
                                const double& mpi,
                                const bool& regularize,
                                const double& R)
  {
    auto integrand =
      [&](double r)
      {
        return (ApplyRegulator(regularize, r, R)
                * WPi(r, mpi) * YPi(r, mpi)
                * std::exp(-square(r / b))
                * std::pow(r / b, l + lp + 2)
                * sf::Laguerre(n, l + 0.5, square(r / b))
                * sf::Laguerre(np, lp + 0.5, square(r / b)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm_product = (ho::CoordinateSpaceNorm(n, l, 1)
                           * ho::CoordinateSpaceNorm(np, lp, 1) / b);
    integral *= norm_product;
    return integral;
  }

  // Integrates m_π r W_π(r) Y_π(r), sandwiched between coordinate space HO
  // radial basis states.
  static double IntegralMPiRWPiRYPiR(const int& np,
                                     const int& lp,
                                     const int& n,
                                     const int& l,
                                     const double& b,
                                     const double& mpi,
                                     const bool& regularize,
                                     const double& R)
  {
    auto integrand =
      [&](double r)
      {
        return (ApplyRegulator(regularize, r, R)
                * std::exp(-square(r / b) - (mpi * r))
                * (3 + 3 * (mpi * r) + square(mpi * r))
                * std::pow(r / b, l + lp)
                * sf::Laguerre(n, l + 0.5, square(r / b))
                * sf::Laguerre(np, lp + 0.5, square(r / b)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm_product = (ho::CoordinateSpaceNorm(n, l, b)
                           * ho::CoordinateSpaceNorm(np, lp, b) / square(mpi));
    integral *= norm_product;
    return integral;
  }

  // Integrates Z_π(r)Y_π(r), sandwiched between coordinate space HO radial
  // basis states.
  static double IntegralZPiYPiR(const int& np,
                                const int& lp,
                                const int& n,
                                const int& l,
                                const double& b,
                                const double& mpi,
                                const bool& regularize,
                                const double& R)
  {
    auto integrand =
      [&](double r)
      {
        return (ApplyRegulator(regularize, r, R)
                * ZPi(r, mpi) * YPi(r, mpi)
                * std::exp(-square(r / b))
                * std::pow(r / b, l + lp + 2)
                * sf::Laguerre(n, l + 0.5, square(r / b))
                * sf::Laguerre(np, lp + 0.5, square(r / b)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm_product = (ho::CoordinateSpaceNorm(n, l, 1)
                           * ho::CoordinateSpaceNorm(np, lp, 1) / b);
    integral *= norm_product;
    return integral;
  }

  // Integrates T_π(r)Y_π(r), sandwiched between coordinate space HO
  // radial basis states.
  static double IntegralTPiYPiR(const int& np,
                                const int& lp,
                                const int& n,
                                const int& l,
                                const double& b,
                                const double& mpi,
                                const bool& regularize,
                                const double& R)
  {
    auto integrand =
      [&](double r)
      {
        return (ApplyRegulator(regularize, r, R)
                * TPi(r, mpi) * YPi(r, mpi)
                * std::exp(-square(r / b))
                * std::pow(r / b, l + lp + 2)
                * sf::Laguerre(n, l + 0.5, square(r / b))
                * sf::Laguerre(np, lp + 0.5, square(r / b)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm_product = (ho::CoordinateSpaceNorm(n, l, 1)
                           * ho::CoordinateSpaceNorm(np, lp, 1) / b);
    integral *= norm_product;
    return integral;
  }

  // Integrates the regularized contact term in the radial momentum space HO basis.
  // Angular momentum selection rules forces l = 0.
  static double IntegralRegularizedDelta(const int& n,
                                         const double& b,
                                         const double& R)
  {
    auto integrand =
      [&](double p)
      {
        return (std::exp(-((square(b * p) / 2) + square(R * p / 2)))
                * square(p) * sf::Laguerre(n, 0.5, square(b * p)));
      };

    double a = 0;
    double epsabs = 0;
    double epsrel = 1e-7;
    double integral = QAGIUIntegrate(integrand, a, epsabs, epsrel);
    double norm = ho::MomentumSpaceNorm(n, 0, b);
    integral *= norm;
    return integral;
  }

  // The regulator function is e^{(p^2 + p'^2) / Λ^2} where Λ is twice the
  // value of the LENPIC regulator.
  // static double F(std::size_t n, double regulator, double result)
  // {
  //   if (n == 0)
  //     return result;
  //   else
  //     return F(n - 1, regulator,
  //              (((2 - square(regulator)) / (2 + square(regulator)))
  //               * std::sqrt((2 * n + 1.0) / (2 * n)) * result));
  // }
  // static double IntegralRegularizedDelta(gsl_params_2n& p)
  // {
  //   auto F_product = (F(p.n, p.scaled_R, 1)
  //                     * F(p.np, p.scaled_R, 1));
  //   auto result = 8 * F_product;
  //   result /= cube(constants::sqrtpi * (2 + p.scaled_R));
  //   return result;
  // }
}

#endif
