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
                                     const double& osc_b)
  {
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////// NLO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b);
  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b);

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                      const basis::RelativeStateLSJT& ket,
                                      const double& osc_b)
  {
    return NLO1Body(bra, ket, osc_b) + NLO2Body(bra, ket, osc_b);
  }

  double NLO1Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b)
  {
    int ni = ket.N(), nf = bra.N();
    int li = ket.L(), lf = bra.L();
    int si = ket.S(), sf = bra.S();
    int ji = ket.J(), jf = bra.J();
    int ti = ket.T(), tf = bra.T();

    // Radial quanta.
    auto nri = (ni - li) / 2;
    auto nrf = (nf - lf) / 2;

    bool kronecker = (nri == nrf && li == lf && si == sf && ti == tf);
    if (!kronecker)
      return 0;

    auto term1 = ((std::sqrt(li * (li + 1.0) * (2 * li + 1.0)) / 2)
                  * std::pow(-1, ji + si + li + 1)
                  * am::Wigner6J(li, jf, si, ji, li, 1));

    auto term2 = (std::sqrt(6)
                  * constants::isoscalar_nucleon_magnetic_moment
                  * std::pow(-1, li + jf)
                  * am::Wigner6J(1, jf, li, ji, 1, 1)
                  * (si == 1 && sf == 1));

    auto result = Hat(ji) * (term1 + term2);
    return  result;
  }

  double NLO2Body(const basis::RelativeStateLSJT& bra,
                  const basis::RelativeStateLSJT& ket,
                  const double& osc_b)
  {
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////////// N2LO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double M1Operator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const double& osc_b)
  {
    return 0;
  }

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////////// N3LO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double IntegralA(int np, int lp, int n, int l, double scale);
  double IntegralB(int np, int lp, int n, int l, double scale);
  double Chi(int n);

  double M1Operator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const double& osc_b)
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

    auto overall_prefactor = (2 * constants::nucleon_mass_fm
                              / std::pow(osc_b, 3));

    // One pion exchange term, with d9.
    auto mpi_b = constants::pion_mass_fm * osc_b;
    auto term1 = (IntegralA(nrf, lf, nri, li, mpi_b)
                  * am::LSCoupledTotalSpinRME(lf, sf, jf, li, si, ji));
    auto term2 = (IntegralB(nrf, lf, nri, li, mpi_b)
                  * am::LSCoupledSDotRhatRhatRME(lf, sf, jf, li, si, ji));
    auto one_pion_exchange_term = term1 - term2;
    one_pion_exchange_term *= (am::PauliDotProductRME(ti, tf)
                               * std::pow(mpi_b, 3) / (4 * constants::pi));
    auto one_pion_exchange_prefactor = (-4 * constants::gA * constants::d9_fm
                                        / std::pow(constants::pion_decay_constant_fm, 2));
    one_pion_exchange_term *= one_pion_exchange_prefactor;

    auto result = overall_prefactor * one_pion_exchange_term;

    // Contact term with d9 and L2.
    bool contact_term_kronecker = (li == 0 && lf == 0 && ji == 1 && jf == 1);
    if (contact_term_kronecker)
      {
        auto d9_contact_term = ((one_pion_exchange_prefactor / 3)
                                * am::PauliDotProductRME(ti, tf));
        auto L2_contact_term = (2 * constants::L2_fm * constants::sqrt2
                                / constants::pi / constants::sqrtpi);
        auto contact_term = ((d9_contact_term + L2_contact_term)
                             * Chi(nrf) * Chi(nri));

        result += overall_prefactor * contact_term;
      }

    return result;
  }

  double Chi(int n)
  {
    return n == 0 ? 1 : std::sqrt((2 * n + 1.0) / (2 * n)) * Chi(n - 1);
  }

  double IntegralA(int np, int lp, int n, int l, double scale)
  {
    if (lp != l)
      return 0;

    struct params { int np_; int n_; int lp_; int l_; double scale_; };
    params p = { np, n, lp, l, scale };
    auto integrand =
      [&p](double y)
      {
        return ((p.lp_ == 0 && p.l_ == 0)
                ? ((std::exp(-p.scale_ * std::sqrt(y)) / y)
                   * gsl_sf_laguerre_n(p.np_, 0.5, y)
                   * gsl_sf_laguerre_n(p.n_, 0.5, y)
                   * (1 - std::exp(-p.scale_ * std::sqrt(y))
                      + p.scale_ * std::sqrt(y)))
                : ((std::exp(-p.scale_ * std::sqrt(y)) / y)
                   * gsl_sf_laguerre_n(p.np_, p.lp_ + 0.5, y)
                   * gsl_sf_laguerre_n(p.n_, p.l_ + 0.5, y)
                   * (1 + p.scale_ * std::sqrt(y))));
      };

    double a = 0;
    double b = 1;
    double alpha = (lp + l) / 2.0;
    int nodes = (np + n) / 2 + 1;
    auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

    auto norm_product = (threedho::CoordinateSpaceNorm(n, l, 1)
                         * threedho::CoordinateSpaceNorm(np, lp, 1));
    norm_product /= (2 * std::pow(scale, 3));

    auto result = norm_product * integral;
    return integral;
  }

  double IntegralB(int np, int lp, int n, int l, double scale)
  {
    struct params { int np_; int n_; int lp_; int l_; double scale_; };
    params p = { np, n, lp, l, scale };
    auto integrand =
      [&p](double y)
      {
        return ((p.lp_ == 0 && p.l_ == 0)
                ? ((std::exp(-p.scale_ * std::sqrt(y)) / y)
                   * gsl_sf_laguerre_n(p.np_, 0.5, y)
                   * gsl_sf_laguerre_n(p.n_, 0.5, y)
                   * (3 - 3 * std::exp(-p.scale_ * std::sqrt(y))
                      + 3 * p.scale_ * std::sqrt(y)
                      + p.scale_ * p.scale_ * y))
                : ((std::exp(-p.scale_ * std::sqrt(y)) / y)
                   * gsl_sf_laguerre_n(p.np_, p.lp_ + 0.5, y)
                   * gsl_sf_laguerre_n(p.n_, p.l_ + 0.5, y)
                   * (3 + 3 * p.scale_ * std::sqrt(y)
                      + p.scale_ * p.scale_ * y)));
      };

    double a = 0;
    double b = 1;
    double alpha = (lp + l) / 2.0;
    int nodes = (np + n) / 2 + 1;
    auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

    auto norm_product = (threedho::CoordinateSpaceNorm(n, l, 1)
                         * threedho::CoordinateSpaceNorm(np, lp, 1));
    norm_product /= (2 * std::pow(scale, 3));

    auto result = norm_product * integral;
    return integral;
  }

  ///////////////////////////////////////////////////////////////////////////////
  /////////////////////////// N4LO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double M1Operator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                       const basis::RelativeStateLSJT& ket,
                                       const double& osc_b)
  {
    return 0;
  }
}
