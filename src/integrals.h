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

  // Integrates the Yukawa function sandwiched between coordinate space radial
  // basis states.
  // static double IntegralYPiR(const gsl_params_2n& p)
  // {
  //   auto integrand =
  //     [&p](double y)
  //     {
  //       return (ApplyRegulator(p.regularize, y, p.scaled_R)
  //               * YPi(y, p.scaled_mpi)
  //               * sf::Laguerre(p.n, p.l + 0.5, y)
  //               * sf::Laguerre(p.np, p.lp + 0.5, y));
  //     };

  //   double a = 0;
  //   double b = 1;
  //   double alpha = (p.l + p.lp + 1.0) / 2.0;
  //   std::size_t nodes = (p.n + p.np) / 2 + 1;
  //   double integral = GaussLaguerre(integrand, a, b, alpha, nodes);
  //   double norm_product = (0.5 * ho::CoordinateSpaceNorm(p.n, p.l, 1)
  //                          * ho::CoordinateSpaceNorm(p.np, p.lp, 1));
  //   integral *= norm_product;
  //   return integral;
  // }

  // // Integrates W_π(r) Y_π(r), sandwiched between coordinate space radial basis
  // // states.
  // static double IntegralWPiRYPiR(const gsl_params_2n& p)
  // {
  //   auto integrand =
  //     [&p](double y)
  //     {
  //       return (ApplyRegulator(p.regularize, y, p.scaled_R)
  //               * WPiYPi(y, p.scaled_mpi)
  //               * sf::Laguerre(p.n, p.l + 0.5, y)
  //               * sf::Laguerre(p.np, p.lp + 0.5, y));
  //     };

  //   double a = 0;
  //   double b = 1;
  //   double alpha = (p.l + p.lp + 1.0) / 2.0;
  //   std::size_t nodes = (p.n + p.np) / 2 + 1;
  //   double integral = GaussLaguerre(integrand, a, b, alpha, nodes);
  //   double norm_product = (0.5 * ho::CoordinateSpaceNorm(p.n, p.l, 1)
  //                          * ho::CoordinateSpaceNorm(p.np, p.lp, 1));
  //   integral *= norm_product;
  //   return integral;
  // }

  // // Integrates m_π r W_π(r) Y_π(r), sandwiched between coordinate space radial
  // // basis states.
  // static double IntegralMPiRWPiRYPiR(const gsl_params_2n& p)
  // {
  //   auto integrand =
  //     [&p](double y)
  //     {
  //      return (ApplyRegulator(p.regularize, y, p.scaled_R)
  //              * WPi(y, p.scaled_mpi)
  //              * std::exp(-std::sqrt(y) * p.scaled_mpi)
  //              * sf::Laguerre(p.n, p.l + 0.5, y)
  //              * sf::Laguerre(p.np, p.lp + 0.5, y));
  //     };

  //   double a = 0;
  //   double b = 1;
  //   double alpha = (p.l + p.lp + 1.0) / 2.0;
  //   std::size_t nodes = (p.n + p.np) / 2 + 1;
  //   auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
  //   return integral;
  // }

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

  // Integrates the regularized contact term in the radial momentum basis.
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
