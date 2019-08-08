#include <cmath>
#include "basis/lsjt_scheme.h"
#include "rme_extras.h"
#include "constants.h"
#include "utility"
#include "integrals.h"
#include "threedho.h"
#include "gt.h"

namespace chiral
{
  // Leading order matrix element.
  double GTOperator::LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator)
  {
    std::size_t nr = ket.n(), nrp = bra.n();
    std::size_t lr = ket.L(), lrp = bra.L();
    std::size_t s = ket.S(), sp = bra.S();
    std::size_t j = ket.J(), jp = bra.J();
    std::size_t t = ket.T(), tp = bra.T();

    bool kronecker = (nr == nrp && lr == lrp);
    if (!kronecker)
      return 0;

    auto symm_term_isospin = 0.5 * am::SpinSymmetricRME(tp, t);
    auto symm_term_spin = (Hat(j) * Hat(s) * ParitySign(jp + lr)
                           * am::Wigner6J(1, j, lr, jp, 1, 1)
                           * am::SpinSymmetricRME(sp, s));
    auto symm_term = symm_term_spin * symm_term_isospin;

    auto asymm_term_isospin = 0.5 * am::SpinAsymmetricRME(tp, t);
    auto asymm_term_spin = (Hat(j) * Hat(sp) * ParitySign(jp + s + lr + 1)
                            * am::Wigner6J(s, j, lr, jp, sp, 1)
                            * am::SpinAsymmetricRME(sp, s));

    auto asymm_term = asymm_term_spin * asymm_term_isospin;

    auto result = (constants::gA * (symm_term + asymm_term));
    if (isnan(result))
      result = 0;
    return result;
  }
}
