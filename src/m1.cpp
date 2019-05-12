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

  double IntegralA(int np, int lp, int n, int l);
  double IntegralB(int np, int lp, int n, int l);
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

    bool kronecker = (si == 1 && sf == 1 && ti == tf);
    if (!kronecker)
      return 0;

    auto prefactor = 2 * constants::nucleon_mass_fm / std::pow(osc_b, 3);

    auto mpi_b = constants::pion_mass_fm * osc_b;
    auto term1 = (IntegralA(nrf, lf, nri, li)
                  * am::LSCoupledTotalSpinRME(lf, sf, jf, li, si, ji));
    auto term2 = (IntegralB(nrf, lf, nri, li)
                  * am::LSCoupledSDotRhatRhatRME(lf, sf, jf, li, si, ji));
    auto term3 = term1 - term2;
    term3 *= (am::PauliDotProductRME(ti, tf)
              * std::pow(mpi_b, 3) / (4 * constants::pi));
    auto term3_prefactor = (-4 * constants::gA * constants::d9_fm
                            / std::pow(constants::pion_decay_constant_fm, 2));
    term3 *= term3_prefactor;

    auto result = prefactor * term3;

    bool term4_kronecker = (li == 0 && lf == 0 && ji == 1 && jf == 1);
    if (term4_kronecker)
      {
        auto term4 = (2 * constants::L2_fm * constants::sqrt2
                      / constants::pi / constants::sqrtpi);
        term4 *= Chi(nrf) * Chi(nri);

        result += prefactor * term4;
      }

    return result;
  }

  double Chi(int n)
  {
    return n == 0 ? 1 : std::sqrt((2 * n + 1.0) / (2 * n)) * Chi(n - 1);
  }

  double IntegralA(int np, int lp, int n, int l)
  {
    return 0;
  }

  double IntegralB(int np, int lp, int n, int l)
  {
    return 0;
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
