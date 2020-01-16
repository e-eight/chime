#ifndef TPRME_H
#define TPRME_H

#include <cstdlib>
#include <cmath>
#include <tuple>
#include "am/wigner_gsl.h"
#include "am/racah_reduction.h"
#include "am/rme.h"
#include "utility.h"

namespace am
{
  namespace tp
  {
    // Container to store am quantum numbers.
    using RelativeState = std::tuple<int, int, int>;
    using RelativeCMState = std::tuple<int, int, int, int, int>;

    // ⟨S'||[σ1 ⊗ σ2]_a||S⟩, a = 0, 1, 2.
    static double SpinProductRME(const int& Sp, const int& S, const int& rank)
    {
      switch (rank)
        {
        case 0:
          return ParitySign(S) * Hat(1 - S) * (1 - std::abs(Sp - S)) / Hat(S);
        case 1:
          return std::sqrt(2) * Hat(S) * std::abs(Sp - S);
        case 2:
          return std::sqrt(20.0 / 3) * (1 - std::abs(Sp - S)) * (1 - std::abs(Sp - 1));
        }
    }

    // ⟨L'S'J'||[C_a(\hat{r}) ⊗ Σ_b]_c||LSJ⟩, where C_a(\hat{r}) are the
    // renormalised spherical harmonics acting on the relative coordinate and
    // Σ_b is a spin operator of rank b.
    template <class SpinOperator>
    static double RelativeSpinRME(const RelativeState& bra,
                                  const RelativeState& ket,
                                  const int& a,
                                  const SpinOperator& Sigma,
                                  const int& b,
                                  const int& c)
    {
      int Lp = std::get<0>(bra), Sp = std::get<1>(bra), Jp = std::get<2>(bra);
      int L = std::get<0>(ket), S = std::get<1>(ket), J = std::get<1>(ket);

      double result = HatProduct(Lp, Sp, J, c);
      result *= am::Wigner9J(L, S, J, a, b, c, Lp, Sp, Jp);
      result *= am::SphericalHarmonicCRME(Lp, L, a);
      result *= Sigma(Sp, S, b);
      return result;
    }

    // ⟨lr' lc' L' S' J'||[[C_a(\hat{r}) ⊗ C_b(\hat{R})]_c ⊗ Σ_d]_e||lr lc L S J⟩.
    // C_b(\hat{R}) are the renormalised spherical harmonics acting on the CM
    // coordinate.
    template <class SpinOperator>
    static double RelativeCMSpinRME(const RelativeCMState& bra,
                                    const RelativeCMState& ket,
                                    const int& a,
                                    const int& b,
                                    const int& c,
                                    const SpinOperator& Sigma,
                                    const int& d,
                                    const int& e)
    {
      int lrp = std::get<0>(bra), lcp = std::get<1>(bra);
      int Lp = std::get<2>(bra), Sp = std::get<3>(bra), Jp = std::get<4>(bra);

      int lr = std::get<0>(ket), lc = std::get<1>(ket);
      int L = std::get<2>(ket), S = std::get<3>(ket), J = std::get<4>(ket);

      double result = HatProduct(Lp, Sp, J, e, lrp, lcp, L, c);
      result *= (am::Wigner9J(L, S, J, c, d, e, Lp, Sp, Jp)
                 * am::Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
      result *= (am::SphericalHarmonicCRME(lrp, lr, a)
                 * am::SphericalHarmonicCRME(lcp, lc, b));
      result *= Sigma(Sp, S, d);
      return result;
    }
  }
}


#endif TPRME_H
