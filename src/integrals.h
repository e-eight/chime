#ifndef INTEGRALS_H
#define INTEGRALS_H

/*
 * Evaluates some of the common integrals that appear in the matrix element
 * evaluation of most operators. The radial integrals are evaluated using
 * Gauss-Laguerre quadrature. For the radial integrals it is assumed that
 * oscillator constant is 1. This requires that parameters like the regulator
 * constant, or the pion mass, be appropriately scaled.
 */

#include "quadrature.h"
#include "constants.h"
#include "specfunc.h"
#include "utility.h"

namespace quadrature
{
  // Structs for passing parameters to GSL integration routines.
  struct gsl_params_2n
  {
    std::size_t np;
    std::size_t lp;
    std::size_t n;
    std::size_t l;
    bool regularize;
    double scaled_regulator;
    double scaled_pion_mass;
  };


  // Integrates m_π R (or m_π r) sandwiched between the radial basis states,
  // in coordinate space. The radial basis states are not normalized.
  static double dd(std::size_t np, std::size_t n, std::size_t lp, std::size_t l)
  {
    return std::sqrt(n + l + 1.5) * (np == n) - std::sqrt(n) * (np + 1 == n);
  }
  static double IntegralMPiR(std::size_t np, std::size_t n, std::size_t lp, std::size_t l)
  {
    double result = 0;
    if (lp == l + 1)
      {
        result = (dd(np, n, lp, l) * am::SphericalHarmonicCRME(lp, l, 1));
      }
    if (l == lp + 1)
      {
        result = (dd(n, np, l, lp) * am::SphericalHarmonicCRME(lp, l, 1));
      }
    return result;
  }

  // Integrates the Yukawa function sandwiched between radial basis states,
  // in coordinate space. The basis states are not normalized.
  static double IntegralYPiR(const gsl_params_2n& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (util::ApplyRegulator(p.regularize,
                                     std::sqrt(y), p.scaled_regulator)
                * std::exp(-p.scaled_pion_mass * std::sqrt(y))
                * sf::Laguerre(p.n, p.l + 0.5, y)
                * sf::Laguerre(p.np, p.lp + 0.5, y));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    auto result = 0.5 * integral / p.scaled_pion_mass;
    return result;
  }

  // Integrates W_π(r) Y_π(r), sandwiched between unnormalized coordinate
  // space radial basis states.
  static double IntegralWPiRYPiR(const gsl_params_2n& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (util::ApplyRegulator(p.regularize,
                                     std::sqrt(y), p.scaled_regulator)
                * util::WPi(std::sqrt(y), p.scaled_pion_mass)
                * std::exp(-p.scaled_pion_mass * std::sqrt(y))
                * sf::Laguerre(p.n, p.l + 0.5, y)
                * sf::Laguerre(p.np, p.lp + 0.5, y));
      };

    double error;
    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    auto result = 0.5 * integral / p.scaled_pion_mass;
    return result;
  }

  // Integrates m_π r W_π(r) Y_π(r), sandwiched between unnormalized coordinate
  // space radial basis states.
  static double IntegralMPiRWPiRYPiR(const gsl_params_2n& p)
  {
    auto integrand =
      [&p](double y)
      {
       return (util::ApplyRegulator(p.regularize,
                                     std::sqrt(y), p.scaled_regulator)
               * util::WPi(std::sqrt(y), p.scaled_pion_mass)
               * std::exp(-p.scaled_pion_mass * std::sqrt(y))
               * sf::Laguerre(p.n, p.l + 0.5, y)
               * sf::Laguerre(p.np, p.lp + 0.5, y));
      };

    /* double error; */
    /* double a = 0; */
    /* double b = 1; */
    /* double alpha = (p.l + p.lp) / 2.0; */
    /* std::size_t nodes = (p.n + p.np) / 2 + 1; */
    /* double integral = gauss_kronrod<double, 15>::integrate(integrand, 1, std::numeric_limits<double>::infinity(), 0, 0, &error); */
    /* return integral; */

    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp + 1.0) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    return integral;
  }

  // Integrates Z_π(r)Y_π(r), sandwiched between unnormalized coordinate space
  // radial basis states.
  static double IntegralZPiYPiR(const gsl_params_2n& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (util::ApplyRegulator(p.regularize,
                                     std::sqrt(y), p.scaled_regulator)
                * util::ZPi(std::sqrt(y), p.scaled_pion_mass)
                * std::exp(-p.scaled_pion_mass * std::sqrt(y))
                * sf::Laguerre(p.n, p.l + 0.5, y)
                * sf::Laguerre(p.np, p.lp + 0.5, y));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    auto result = 0.5 * integral / p.scaled_pion_mass;
    return result;
  }

  // Integrates T_π(r)Y_π(r), sandwiched between unnormalized coordinate space
  // radial basis states.
  static double IntegralTPiYPiR(const gsl_params_2n& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (util::ApplyRegulator(p.regularize,
                                     std::sqrt(y), p.scaled_regulator)
                * util::TPi(std::sqrt(y), p.scaled_pion_mass)
                * std::exp(-p.scaled_pion_mass * std::sqrt(y))
                * sf::Laguerre(p.n, p.l + 0.5, y)
                * sf::Laguerre(p.np, p.lp + 0.5, y));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    auto result = 0.5 * integral / p.scaled_pion_mass;
    return result;
  }

  // Integrates the regularized contact term in the radial momentum basis.
  // The regulator function is e^{(p^2 + p'^2) / Λ^2} where Λ is twice the
  // value of the LENPIC regulator.
  static double F(std::size_t n, double regulator, double result)
  {
    if (n == 0)
      return result;
    else
      return F(n - 1, regulator,
               (((2 - square(regulator)) / (2 + square(regulator)))
                * std::sqrt((2 * n + 1.0) / (2 * n)) * result));
  }
  static double IntegralRegularizedDelta(gsl_params_2n& p)
  {
    auto F_product = (F(p.n, p.scaled_regulator, 1)
                      * F(p.np, p.scaled_regulator, 1));
    auto result = 8 * F_product;
    result /= cube(constants::sqrtpi * (2 + p.scaled_regulator));
    return result;
  }
}

#endif
