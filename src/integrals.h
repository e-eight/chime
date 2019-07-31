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
  // Struct for passing parameters to GSL integration routines.
  struct gsl_params
  {
    std::size_t n;
    std::size_t l;
    std::size_t np;
    std::size_t lp;
    bool regularize;
    double scaled_regulator;
    double scaled_pion_mass;
  };

  // Integrates m_π R (or m_π r) sandwiched between the radial basis states,
  // in coordinate space. The radial basis states are not normalized.
  double IntegralMPiR(const gsl_params& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (util::ApplyRegulator(p.regularize,
                                     std::sqrt(y), p.scaled_regulator)
                * sf::Laguerre(p.n, p.l + 0.5, y)
                * sf::Laguerre(p.np, p.lp + 0.5, y));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp + 2.0) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    auto result = 0.5 * p.scaled_pion_mass * integral;
    return result;
  }

  // Integrates the Yukawa function sandwiched between radial basis states,
  // in coordinate space. The basis states are not normalized.
  double IntegralYPiR(const gsl_params& p)
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

  // Integrates m_π r W_π(r) Y_π(r), sandwiched between unnormalized coordinate
  // space radial basis states.
  double IntegralMPiRWPiRYPiR(const gsl_params& p)
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

    double a = 0;
    double b = 1;
    double alpha = (p.l + p.lp + 1.0) / 2.0;
    std::size_t nodes = (p.n + p.np) / 2 + 1;
    auto integral = GaussLaguerre(integrand, a, b, alpha, nodes);
    return integral;
  }

  // Integrates Z_π(r)Y_π(r), sandwiched between unnormalized coordinate space
  // radial basis states.
  double IntegralZPiYPiR(const gsl_params& p)
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
  double IntegralTPiYPiR(const gsl_params& p)
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

}

#endif
