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
#include "m1.h"

using namespace constants;
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

  double M1Operator::NLOMatrixElement(const basis::RelativeStateLSJT& bra,
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

  double IntegralA(int ni, int li, int nf, int lf);
  double IntegralB(int ni, int li, int nf, int lf);
  double SpaceSpinDotRME(int li, int si, int ji, int lf, int sf, int jf);
  double Chi(int n, double result);

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

    auto mpi_b = constants::pion_mass_fm * osc_b;
    auto term1 = (IntegralA(ni, li, nf, lf)
                  * util::TotalSpin(li, si, ji, lf, sf, jf));
    auto term2 = (IntegralB(ni, li, nf, lf)
                  * SpaceSpinDotRME(li, si, ji, lf, sf, jf));
    auto term3 = term1 - term2;
    term3 *= (util::PauliDotProduct(ti, tf)
              * std::pow(mpi_b, 3) / (4 * constants::pi));
    auto term3_prefactor = (-4 * constants::gA * constants::d9_fm
                            / std::pow(constants::pion_decay_constant_fm, 2));
    term3 *= term3_prefactor;

    auto term4 = (2 * constants::L2_fm * constants::sqrt2
                  / constants::pi / constants::sqrtpi);
    term4 *= (li == 0 && lf == 0 && ji == 1 && jf == 1);
    term4 *= Chi(nf, 1) * Chi(ni, 1);

    auto prefactor = 2 * constants::nucleon_mass_fm / std::pow(osc_b, 3);

    auto result = prefactor * (term3 + term4);
    return result;
  }

  double Chi(int n, double result)
  {
    if (n == 1)
      return result;
    return Chi(n - 1, std::sqrt((2 * n + 1.0) / n) * result);
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
