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
