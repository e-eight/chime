#ifndef RME_EXTRAS_H
#define RME_EXTRAS_H

/*
 * Defines reduced matrix elements for two-body Pauli matrices,
 * and their tensor products with unnormalized spherical harmonics.
 */


#include <cmath>
#include "am/wigner_gsl.h"
#include "am/racah_reduction.h"
#include "am/rme.h"

namespace am
{
  // Calculate reduced matrix element of Pauli matrices in two-body (iso)spin
  // space.
  //
  // Arguments:
  //  sp (int): bra spin
  //  s (int): ket spin
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double PauliOneRME(const int& sp, const int& s)
  {
    assert(AllowedTriangle(HalfInt(1, 2), HalfInt(1, 2), sp));
    assert(AllowedTriangle(HalfInt(1, 2), HalfInt(1, 2), s));

    double result = 0;
    if (sp == s)
      result = std::sqrt(s * (s + 1.));
    else
      result = ParitySign(s) * std::sqrt(3);
    return result;
  }
  inline double PauliTwoRME(const int& sp, const int& s)
  {
    assert(AllowedTriangle(HalfInt(1, 2), HalfInt(1, 2), sp));
    assert(AllowedTriangle(HalfInt(1, 2), HalfInt(1, 2), s));

    double result = 0;
    if (sp == s)
      result = std::sqrt(s * (s + 1.));
    else
      result = ParitySign(sp) * std::sqrt(3);
    return result;
  }

  // Calculate reduced matrix element of [\vec{σ}_1 ⊗ \vec{σ}_2]_k in a
  // two-body coupled (iso)spin basis.
  //
  // Arguments:
  //  sp (int): bra spin
  //  s (int): ket spin
  //  k (int): rank of the tensor product;
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double PauliProductRME(const int& sp, const int& s, const int& k)
  {
    assert(AllowedTriangle(1, 1, k));

    auto hat_product = 2 * Hat(s) * Hat(k);
    auto half = HalfInt(1, 2);
    auto wigner_9j = Wigner9J(half, half, s, 1, 1, k, half, half, sp);
    auto result = 3 * hat_product * wigner_9j;

    return result;
  }

  // Calculate reduced matrix element of
  // [C_a(\hat{r}) ⊗ [\vec{σ}_1 ⊗ \vec{σ}_2]_b]_c in a LS coupled basis.
  //
  // Arguments:
  //  lp (int): bra orbital angular momentum
  //  l (int): ket orbital angular momentum
  //  sp (int): bra spin
  //  s (int): ket spin
  //  jp (int): bra angular momentum
  //  j (int): ket angular momentum
  //  a (int): rank of C_a
  //  b (int): rank of Pauli tensor product
  //  c (int): rank of tensor product
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double RelativePauliProductRME(const int& lp,
                                        const int& l,
                                        const int& sp,
                                        const int& s,
                                        const int& jp,
                                        const int& j,
                                        const int& a,
                                        const int& b,
                                        const int& c)
  {
    assert(AllowedTriangle(a, b, c));

    auto hat_product = Hat(lp) * Hat(sp) * Hat(j) * Hat(c);
    auto wigner_9j = Wigner9J(l, s, j, a, b, c, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto pprme = PauliProductRME(sp, s, b);
    auto result = (hat_product * wigner_9j * crme * pprme);
    return result;
  }

  // Calculate reduced matrix element of
  // [C_a(\hat{r}) ⊗ \vec{σ}_i]_b in a LS coupled basis.
  //
  // Arguments:
  //  lp (int): bra orbital angular momentum
  //  l (int): ket orbital angular momentum
  //  sp (int): bra spin
  //  s (int): ket spin
  //  jp (int): bra angular momentum
  //  j (int): ket angular momentum
  //  a (int): rank of C_a
  //  b (int): rank of tensor product
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double RelativePauliOneRME(const int& lp,
                                    const int& l,
                                    const int& sp,
                                    const int& s,
                                    const int& jp,
                                    const int& j,
                                    const int& a,
                                    const int& b)
  {
    assert(AllowedTriangle(a, 1, b));

    auto hat_product = Hat(lp) * Hat(sp) * Hat(j) * Hat(b);
    auto wigner_9j = Wigner9J(l, s, j, a, 1, b, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto prme = PauliOneRME(sp, s);
    auto result = (hat_product * wigner_9j * crme * prme);
    return result;
  }
  inline double RelativePauliTwoRME(const int& lp,
                                    const int& l,
                                    const int& sp,
                                    const int& s,
                                    const int& jp,
                                    const int& j,
                                    const int& a,
                                    const int& b)
  {
    assert(AllowedTriangle(a, 1, b));

    auto hat_product = Hat(lp) * Hat(sp) * Hat(j) * Hat(b);
    auto wigner_9j = Wigner9J(l, s, j, a, 1, b, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto prme = PauliTwoRME(sp, s);
    auto result = (hat_product * wigner_9j * crme * prme);
    return result;
  }
  inline double RelativeTotalSpinRME(const int& lp,
                                     const int& l,
                                     const int& sp,
                                     const int& s,
                                     const int& jp,
                                     const int& j,
                                     const int& a,
                                     const int& b)
  {
    assert(AllowedTriangle(a, 1, b));

    auto hat_product = Hat(lp) * Hat(sp) * Hat(j) * Hat(b);
    auto wigner_9j = Wigner9J(l, s, j, a, 1, b, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto prme = (PauliOneRME(sp, s) + PauliTwoRME(sp, s))/2;
    auto result = (hat_product * wigner_9j * crme * prme);
    return result;
  }

  // Calculate reduced matrix element of
  // [[C_a(\hat{r}) ⊗ C_b(\hat{R})]_c ⊗ [\vec{σ}_1 ⊗ \vec{σ}_2]_d]_e
  // in a LS coupled basis.
  //
  // Arguments:
  //  lrp (int): bra orbital angular momentum, relative
  //  lr (int): ket orbital angular momentum, relative
  //  lcp (int): bra orbital angular momentum, cm
  //  lc (int): ket orbital angular momentum, cm
  //  Lp (int): bra total orbital angular momentum
  //  L (int): ket total orbital angular momentum
  //  Sp (int): bra spin
  //  S (int): ket spin
  //  Jp (int): bra angular momentum
  //  J (int): ket angular momentum
  //  a (int): rank of C_a
  //  b (int): rank of C_b
  //  c (int): rank of C_a C_b tensor product
  //  d (int): rank of Pauli tensor product
  //  e (int): rank of tensor product
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double RelativeCMPauliProductRME(const int& lrp,
                                          const int& lr,
                                          const int& lcp,
                                          const int& lc,
                                          const int& Lp,
                                          const int& L,
                                          const int& Sp,
                                          const int& S,
                                          const int& Jp,
                                          const int& J,
                                          const int& a,
                                          const int& b,
                                          const int& c,
                                          const int& d,
                                          const int& e)
  {
    assert(AllowedTriangle(a, b, c));
    assert(AllowedTriangle(c, d, e));

    auto hat_product = (Hat(Lp) * Hat(Sp) * Hat(J) * Hat(lrp)
                        * Hat(lcp) * Hat(c) * Hat(L));
    auto wigner_product = (Wigner9J(L, S, J, c, d, e, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto pprme = PauliProductRME(Sp, S, d);
    auto result = (hat_product * wigner_product * crme_product * pprme);
    return result;
  }

  // Calculate reduced matrix element of
  // [[C_a(\hat{r}) ⊗ C_b(\hat{R})]_c ⊗ \vec{σ}_i]_d
  // in a LS coupled basis.
  //
  // Arguments:
  //  lrp (int): bra orbital angular momentum, relative
  //  lr (int): ket orbital angular momentum, relative
  //  lcp (int): bra orbital angular momentum, cm
  //  lc (int): ket orbital angular momentum, cm
  //  Lp (int): bra total orbital angular momentum
  //  L (int): ket total orbital angular momentum
  //  Sp (int): bra spin
  //  S (int): ket spin
  //  Jp (int): bra angular momentum
  //  J (int): ket angular momentum
  //  a (int): rank of C_a
  //  b (int): rank of C_b
  //  c (int): rank of C_a C_b tensor product
  //  d (int): rank of tensor product
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double RelativeCMPauliOneRME(const int& lrp,
                                      const int& lr,
                                      const int& lcp,
                                      const int& lc,
                                      const int& Lp,
                                      const int& L,
                                      const int& Sp,
                                      const int& S,
                                      const int& Jp,
                                      const int& J,
                                      const int& a,
                                      const int& b,
                                      const int& c,
                                      const int& d)
  {
    assert(AllowedTriangle(a, b, c));
    assert(AllowedTriangle(c, 1, d));

    auto hat_product = (Hat(Lp) * Hat(Sp) * Hat(J) * Hat(lrp)
                        * Hat(lcp) * Hat(c) * Hat(L));
    auto wigner_product = (Wigner9J(L, S, J, c, 1, d, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto prme = PauliOneRME(Sp, S);
    auto result = (hat_product * wigner_product * crme_product * prme);
    return result;
  }
  inline double RelativeCMPauliTwoRME(const int& lrp,
                                      const int& lr,
                                      const int& lcp,
                                      const int& lc,
                                      const int& Lp,
                                      const int& L,
                                      const int& Sp,
                                      const int& S,
                                      const int& Jp,
                                      const int& J,
                                      const int& a,
                                      const int& b,
                                      const int& c,
                                      const int& d)
  {
    assert(AllowedTriangle(a, b, c));
    assert(AllowedTriangle(c, 1, d));

    auto hat_product = (Hat(Lp) * Hat(Sp) * Hat(J) * Hat(lrp)
                        * Hat(lcp) * Hat(c) * Hat(L));
    auto wigner_product = (Wigner9J(L, S, J, c, 1, d, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto prme = PauliTwoRME(Sp, S);
    auto result = (hat_product * wigner_product * crme_product * prme);
    return result;
  }
  inline double RelativeCMTotalSpinRME(const int& lrp,
                                       const int& lr,
                                       const int& lcp,
                                       const int& lc,
                                       const int& Lp,
                                       const int& L,
                                       const int& Sp,
                                       const int& S,
                                       const int& Jp,
                                       const int& J,
                                       const int& a,
                                       const int& b,
                                       const int& c,
                                       const int& d)
  {
    assert(AllowedTriangle(a, b, c));
    assert(AllowedTriangle(c, 1, d));

    auto hat_product = (Hat(Lp) * Hat(Sp) * Hat(J) * Hat(lrp)
                        * Hat(lcp) * Hat(c) * Hat(L));
    auto wigner_product = (Wigner9J(L, S, J, c, 1, d, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto prme = (PauliOneRME(Sp, S) + PauliTwoRME(Sp, S)) / 2;
    auto result = (hat_product * wigner_product * crme_product * prme);
    return result;
  }

}

#endif
