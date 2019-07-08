#include <cmath>
#include <unordered_map>
#include <gsl/gsl_sf_laguerre.h>
#include "basis/lsjt_scheme.h"
#include "rme_extended.h"
#include "constants.h"
#include "utility.h"
#include "chiral.h"
#include "quadrature.h"
#include "threedho.h"
#include "m1.h"

namespace chiral
{
  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////// LO Matrix Element ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double M1Operator::LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator)
  {
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////// NLO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b,
                  const double& regulator);
  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b,
                  const double& regulator);

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const double& osc_b,
                                      const double& regulator)
  {
    return (NLO1Body(bra, ket, osc_b, regulator)
            + NLO2Body(bra, ket, osc_b, regulator));
  }

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b,
                  const double& regulator)
  {
    int ni = ket.N(), nf = bra.N();
    int li = ket.L(), lf = bra.L();
    int si = ket.S(), sf = bra.S();
    int ji = ket.J(), jf = bra.J();
    int ti = ket.T(), tf = bra.T();

    // Radial quanta.
    auto nri = (ni - li) / 2;
    auto nrf = (nf - lf) / 2;

    bool kronecker = (nri == nrf && li == lf);
    if (!kronecker)
      return 0;

    auto symm_term_angular = ((1 + 0.5 * std::sqrt(ti * (ti + 1)))
                                  * std::pow(-1, 1 + li + si + ji)
                                  * Hat(li) * std::sqrt(li * (li + 1))
                                  * am::Wigner6J(li, ji, si, jf, li, 1));

    auto symm_term_spin = ((constants::isoscalar_nucleon_magnetic_moment
                                + (std::sqrt(ti * (ti + 1))
                                   * constants::isovector_nucleon_magnetic_moment))
                               * std::pow(-1, li + jf)
                               * Hat(si) * std::sqrt(si * (si + 1))
                               * am::Wigner6J(1, ji, li, jf, 1, 1));

    auto symm_term = (symm_term_angular + symm_term_spin) * (tf == ti) * (sf == si);

    auto asymm_term = (0.75 * (ParitySign(ti) - ParitySign(tf))
                       * (ParitySign(si) - ParitySign(sf))
                       * ParitySign(li + si + jf + 1) * Hat(sf)
                       * am::Wigner6J(si, ji, li, jf, sf, 1)
                       * (tf != ti) * (sf != si));

    auto result = Hat(ji) * (symm_term + asymm_term);
    return result;
  }

  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b,
                  const double& regulator)
  {
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////////// N2LO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double M1Operator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const double& osc_b,
                                       const double& regulator)
  {
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////////// N3LO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  struct gsl_params { int nrf; int nri; int lf; int li; double scale; double regulator; };
  double IntegralA(const gsl_params& p);
  double IntegralB(const gsl_params& p);
  double IntegralC(const gsl_params& p);

  double M1Operator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const double& osc_b,
                                       const double& regulator)
  {
    int ni = ket.N(), nf = bra.N();
    int li = ket.L(), lf = bra.L();
    int si = ket.S(), sf = bra.S();
    int ji = ket.J(), jf = bra.J();
    int ti = ket.T(), tf = bra.T();

    // Radial quanta.
    auto nri = (ni - li) / 2;
    auto nrf = (nf - lf) / 2;

    // Overall selection criteria.
    bool kronecker = (si == 1 && sf == 1 && ti == tf);
    if (!kronecker)
      return 0;

    auto mpi_b = constants::pion_mass_fm * osc_b;
    auto regulator_b = regulator / osc_b;
    auto d9_factor = ((4 * constants::d9_fm * constants::gA
                       / std::pow(constants::pion_decay_constant_fm, 2))
                      * am::PauliDotProductRME(tf, ti));
    gsl_params p {nrf, nri, lf, li, mpi_b, regulator_b};

    // One pion exchange term with d9.
    auto term1 = ((std::sqrt(1.0 / 3.0)
                   * am::LSCoupledTotalSpinY0Rank1RME(lf, sf, jf, li, si, ji))
                  - (std::sqrt(2.0 / 3.0)
                     * am::LSCoupledTotalSpinY2Rank1RME(lf, sf, jf, li, si, ji)));
    term1 *= (IntegralA(p) * std::sqrt(4.0 * constants::pi / 3.0));
    auto term2 = IntegralB(p) * am::LSCoupledTotalSpinRME(lf, sf, jf, li, si, ji);
    auto one_pion_exchange_term = d9_factor * (term1 - term2);
    one_pion_exchange_term *= (std::pow(constants::pion_mass_fm, 3)
                               / (4 * constants::pi));

    // Contact term with d9 and L2.
    double contact_term = 0;
    if (li == 0 && lf == 0)
      {
        auto L2_d9_factor = (2 * constants::L2_fm - (d9_factor / 3));
        auto scaled_IntegralC = IntegralC(p);
        scaled_IntegralC /= (4 * std::pow(constants::pi, 1.5)
                             * regulator * osc_b * osc_b);
        contact_term = (L2_d9_factor * scaled_IntegralC
                        * am::LSCoupledTotalSpinRME(lf, sf, jf, li, si, ji));
      }

    // Total value in terms of nuclear magneton.
    auto result = one_pion_exchange_term + contact_term;
    result /= constants::nuclear_magneton_fm;
    return result;
  }

  // double Chi(int n)
  // {
  //   return n == 0 ? 1 : std::sqrt((2 * n + 1.0) / (2 * n)) * Chi(n - 1);
  // }

  double IntegralA(const gsl_params& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (std::exp(-p.scale * std::sqrt(y))
                * gsl_sf_laguerre_n(p.nrf, p.lf + 0.5, y)
                * gsl_sf_laguerre_n(p.nri, p.li + 0.5, y)
                * util::TPi(y, p.scale)
                * util::LenpicSemiLocalRegulator(std::sqrt(y), p.regulator));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.lf + p.li) / 2.0;
    int nodes = (p.nrf + p.nri) / 2 + 1;
    auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

    auto norm_product = (threedho::CoordinateSpaceNorm(p.nri, p.li, 1)
                         * threedho::CoordinateSpaceNorm(p.nrf, p.lf, 1));
    norm_product /= (2 * p.scale);

    auto result = norm_product * integral;
    return integral;
  }

  double IntegralB(const gsl_params& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (std::exp(-p.scale * std::sqrt(y))
                * gsl_sf_laguerre_n(p.nrf, p.lf + 0.5, y)
                * gsl_sf_laguerre_n(p.nri, p.li + 0.5, y)
                * util::ZPi(y, p.scale)
                * util::LenpicSemiLocalRegulator(std::sqrt(y), p.regulator));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.lf + p.li) / 2.0;
    int nodes = (p.nrf + p.nri) / 2 + 1;
    auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

    auto norm_product = (threedho::CoordinateSpaceNorm(p.nri, p.li, 1)
                         * threedho::CoordinateSpaceNorm(p.nrf, p.lf, 1));
    norm_product /= (2 * p.scale);

    auto result = norm_product * integral;
    return integral;
  }

  double IntegralC(const gsl_params& p)
  {
    auto integrand =
      [&p](double y)
      {
        return (std::exp(-y / std::pow(p.regulator, 2))
                * gsl_sf_laguerre_n(p.nrf, p.lf + 0.5, y)
                * gsl_sf_laguerre_n(p.nri, p.li + 0.5, y));
      };

    double a = 0;
    double b = 1;
    double alpha = (p.lf + p.li - 1) / 2.0;
    int nodes = (p.nrf + p.nri) / 2 + 1;
    auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

    auto norm_product = (threedho::CoordinateSpaceNorm(p.nri, p.li, 1)
                         * threedho::CoordinateSpaceNorm(p.nrf, p.lf, 1));

    auto result = norm_product * integral;
    return integral;
  }

  // double IntegralB(int np, int lp, int n, int l, double scale)
  // {
  //   struct params { int np_; int n_; int lp_; int l_; double scale_; };
  //   params p = { np, n, lp, l, scale };
  //   auto integrand =
  //     [&p](double y)
  //     {
  //       // return ((p.lp_ == 0 && p.l_ == 0)
  //       //         ? ((std::exp(-p.scale_ * std::sqrt(y)) / y)
  //       //            * gsl_sf_laguerre_n(p.np_, 0.5, y)
  //       //            * gsl_sf_laguerre_n(p.n_, 0.5, y)
  //       //            * (3 - 3 * std::exp(-p.scale_ * std::sqrt(y))
  //       //               + 3 * p.scale_ * std::sqrt(y)
  //       //               + p.scale_ * p.scale_ * y))
  //       //         : ((std::exp(-p.scale_ * std::sqrt(y)) / y)
  //       //            * gsl_sf_laguerre_n(p.np_, p.lp_ + 0.5, y)
  //       //            * gsl_sf_laguerre_n(p.n_, p.l_ + 0.5, y)
  //       //            * (3 + 3 * p.scale_ * std::sqrt(y)
  //       //               + p.scale_ * p.scale_ * y)));
  //       return ((std::exp(-p.scale_ * std::sqrt(y)) / y)
  //               * gsl_sf_laguerre_n(p.np_, p.lp_ + 0.5, y)
  //               * gsl_sf_laguerre_n(p.n_, p.l_ + 0.5, y)
  //               * (3 + 3 * p.scale_ * std::sqrt(y) + p.scale_ * p.scale_ * y)
  //               * LenpicSemiLocalRegulator(std::sqrt(y), 0.9
  //                                          * constants::pion_mass_fm / p.scale_));
  //     };

  //   double a = 0;
  //   double b = 1;
  //   double alpha = (lp + l) / 2.0;
  //   int nodes = (np + n) / 2 + 1;
  //   auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

  //   auto norm_product = (threedho::CoordinateSpaceNorm(n, l, 1)
  //                        * threedho::CoordinateSpaceNorm(np, lp, 1));
  //   norm_product /= (2 * std::pow(scale, 3));

  //   auto result = norm_product * integral;
  //   return integral;
  // }

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////////// N4LO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double M1Operator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const double& osc_b,
                                       const double& regulator)
  {
    return 0;
  }
}
