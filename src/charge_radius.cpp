#include <cmath>
#include <gsl/gsl_sf_laguerre.h>
#include "basis/lsjt_scheme.h"
#include "basis/am/halfint.h"
#include "basis/am/wigner_gsl.h"
#include "constants.h"
#include "utility.h"
#include "chiral.h"
#include "quadrature.h"
#include "threedho.h"
#include "charge_radius.h"

using namespace constants;
namespace chiral
{
  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////// LO Matrix Element ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double ChargeRadiusOperator::LOMatrixElement(const basis::RelativeStateLSJT& bra,
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

    bool diagonal = (li == lf && si == sf && ji == jf && ti == tf);
    if (!diagonal)
      return 0;

    auto term1 = (2 * nri + li + 1.5) * util::KroneckerDelta(nri, nrf);
    auto term2 = std::sqrt((nri + 1) * (nri + li + 1.5)) * util::KroneckerDelta(nri + 1, nrf);
    auto term3 = std::sqrt((nrf + 1) * (nrf + lf + 1.5)) * util::KroneckerDelta(nrf + 1, nri);
    auto radial_integral = osc_b * osc_b * (term1 - term2 - term3);

    double clebsch_product = 0;
    for (int ms = -si; ms <= si; ++ms)
        clebsch_product += std::pow(am::ClebschGordan(li, -ms, si, ms, ji, 0), 2);

    auto wigner_factor = 1.0;

    auto result = wigner_factor * clebsch_product * radial_integral;
    return result;
  }

  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////// NLO Matrix Element ///////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double ChargeRadiusOperator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
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

    bool diagonal = (nri == nrf && li == lf && si == sf && ji == jf && ti == tf);
    if (!diagonal)
      return 0;

    double clebsch_product = 0;
    for (int ms = -si; ms <= si; ++ms)
        clebsch_product += std::pow(am::ClebschGordan(li, -ms, si, ms, ji, 0), 2);

    auto wigner_factor = 1.0;

    auto result = (wigner_factor * clebsch_product
                   * constants::isoscalar_nucleon_charge_radius_sq_fm);
    return result;
  }


  ///////////////////////////////////////////////////////////////////////////////
  //////////////////////////// N2LO Matrix Element //////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double ChargeRadiusOperator::N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                                 const basis::RelativeStateLSJT& ket,
                                                 const double& osc_b)
  {
    return 0;
  }

    ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////// N3LO Matrix Element ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double ChargeRadiusOperator::N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
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

    bool diagonal = (li == lf && si == sf && ji == jf && ti == tf);
    if (!diagonal)
      return 0;

    // Radial integral
    struct params { int nri_; int nrf_; int li_; int lf_; double mb_; };
    auto mpi_b = constants::pion_mass_fm * osc_b;
    params p = { nri, nrf, li, lf, mpi_b };
    auto integrand =
      [&p](double y) { return (std::exp(-p.mb_ * std::sqrt(y))
                               * gsl_sf_laguerre_n(p.nri_, p.li_ + 0.5, y)
                               * gsl_sf_laguerre_n(p.nrf_, p.lf_ + 0.5, y)); };
    double a = 0;
    double b = 1;
    double alpha = li;
    int nodes = (nri + nrf) / 2 + 1;
    auto integral = quadrature::GaussLaguerre(integrand, a, b, alpha, nodes);

    // Radial normalization
    auto norm_product = (threedho::CoordinateSpaceNorm(nri, li, osc_b)
                           * threedho::CoordinateSpaceNorm(nrf, lf, osc_b));
    norm_product *= (osc_b * osc_b);

    // Spin-Isospin matrix element
    auto spin_element = util::PauliDotProduct(si, sf);
    auto isospin_element = util::PauliDotProduct(ti, tf);

    // Clebsch product
    double clebsch_product = 0;
    for (int ms = -si; ms <= si; ++ms)
      clebsch_product += std::pow(am::ClebschGordan(li, -ms, si, ms, ji, 0), 2);

    // Wigner factor
    auto wigner_factor = 1.0;

    // LEC prefactor
    auto prefactor = -std::pow(constants::gA / constants::pion_decay_constant_fm, 2);
    prefactor /= (64 * constants::pi * constants::nucleon_mass_fm);

    auto result = integral * norm_product * spin_element * isospin_element;
    result *= (prefactor * wigner_factor * clebsch_product);
    return result;
  }

  ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////// N4LO Matrix Element ////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  double ChargeRadiusOperator::N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                                 const basis::RelativeStateLSJT& ket,
                                                 const double& osc_b)
  {
    return 0;
  }
}
