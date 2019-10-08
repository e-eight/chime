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
#include "integrals.h"
#include "threedho.h"

// Product of Hats.
template <class T>
double HatProduct(T x)
{
  return Hat(x);
}
template <class T, class... Args>
double HatProduct(T x, Args... args)
{
  return Hat(x) * HatProduct(args...);
}

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
      result = std::sqrt(s * (s + 1));
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
      result = std::sqrt(s * (s + 1));
    else
      result = ParitySign(sp) * std::sqrt(3);
    return result;
  }

  // Calculate reduced matrix element of (\vec{σ}_1 + \vec{σ}_2) / 2 in its
  // eigenbasis.
  // Arguments:
  //  sp (int): bra spin
  //  s (int): ket spin
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double SpinSymmetricRME(const int& sp, const int& s)
  {
    return AngularMomentumJRME(sp, s);
  }

  // Calculate reduced matrix element of (\vec{σ}_1 - \vec{σ}_2) / 2 in its
  // eigenbasis.
  // Arguments:
  //  sp (int): bra spin
  //  s (int): ket spin
  // Returns:
  //  reduced matrix element (double), Rose convention
  inline double SpinAntisymmetricRME(const int& sp, const int& s)
  {
    return ((s == 0 && sp == 1) - std::sqrt(3) * (s == 1 && sp == 0));
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

    auto hat_product = 2 * HatProduct(s, k);
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

    auto hat_product = HatProduct(lp, sp, j, c);
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

    auto hat_product = HatProduct(lp, sp, j, b);
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

    auto hat_product = HatProduct(lp, sp, j, b);
    auto wigner_9j = Wigner9J(l, s, j, a, 1, b, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto prme = PauliTwoRME(sp, s);
    auto result = (hat_product * wigner_9j * crme * prme);
    return result;
  }
  inline double RelativeSpinSymmetricRME(const int& lp,
                                         const int& l,
                                         const int& sp,
                                         const int& s,
                                         const int& jp,
                                         const int& j,
                                         const int& a,
                                         const int& b)
  {
    // assert(AllowedTriangle(a, 1, b));

    auto hat_product = HatProduct(lp, sp, j, b);
    auto wigner_9j = Wigner9J(l, s, j, a, 1, b, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto prme = SpinSymmetricRME(sp, s);
    auto result = (hat_product * wigner_9j * crme * prme);
    return result;
  }
  inline double RelativeSpinAntisymmetricRME(const int& lp,
                                             const int& l,
                                             const int& sp,
                                             const int& s,
                                             const int& jp,
                                             const int& j,
                                             const int& a,
                                             const int& b)
  {
    assert(AllowedTriangle(a, 1, b));

    auto hat_product = HatProduct(lp, sp, j, b);
    auto wigner_9j = Wigner9J(l, s, j, a, 1, b, lp, sp, jp);
    auto crme = SphericalHarmonicCRME(lp, l, a);
    auto prme = SpinAntisymmetricRME(sp, s);
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

    auto hat_product = HatProduct(Lp, Sp, J, lrp, lcp, c, L);
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

    auto hat_product = HatProduct(Lp, Sp, J, lrp, lcp, c, L);
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

    auto hat_product = HatProduct(Lp, Sp, J, lrp, lcp, c, L);
    auto wigner_product = (Wigner9J(L, S, J, c, 1, d, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto prme = PauliTwoRME(Sp, S);
    auto result = (hat_product * wigner_product * crme_product * prme);
    return result;
  }
  inline double RelativeCMSpinSymmetricRME(const int& lrp,
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

    auto hat_product = HatProduct(Lp, Sp, J, lrp, lcp, c, L);
    auto wigner_product = (Wigner9J(L, S, J, c, 1, d, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto prme = SpinSymmetricRME(Sp, S);
    auto result = (hat_product * wigner_product * crme_product * prme);
    return result;
  }
  inline double RelativeCMSpinAntisymmetricRME(const int& lrp,
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

    auto hat_product = HatProduct(Lp, Sp, J, lrp, lcp, c, L);
    auto wigner_product = (Wigner9J(L, S, J, c, 1, d, Lp, Sp, Jp)
                           * Wigner9J(lr, lc, L, a, b, c, lrp, lcp, Lp));
    auto crme_product = (SphericalHarmonicCRME(lrp, lr, a)
                         * SphericalHarmonicCRME(lcp, lc, b));
    auto prme = SpinAntisymmetricRME(Sp, S);
    auto result = (hat_product * wigner_product * crme_product * prme);
    return result;
  }

  // Calculate reduced matrix element of L_{rel} in a LS coupled basis.
  // Arguments:
  //  lp (int): bra orbital angular momentum, relative
  //  l  (int): ket orbital angular momentum, relative
  //  sp  (int): bra spin
  //  s   (int): ket spin
  //  jp  (int): bra angular momentum
  //  j   (int): ket angular momentum
  // Returns:
  //  reduced matrix element (double)
  inline double RelativeLrelRME(const int& lp, const int& l,
                                const int& sp, const int& s,
                                const int& jp, const int& j)
  {
    if (sp != s)
      return 0;
    auto hat_product = HatProduct(j, lp);
    auto wigner_6j = am::Wigner6J(l, lp, 1, jp, j, s);
    auto parity = ParitySign(lp + s + j + 1);
    auto lrel_rme = AngularMomentumJRME(lp, l);
    auto result = (parity * hat_product * wigner_6j * lrel_rme);
    return result;
  }

  // Calculate reduced matrix element of L_{cm} + L_{rel} in a LS coupled basis.
  // Arguments:
  //  lrp (int): bra orbital angular momentum, relative
  //  lr  (int): ket orbital angular momentum, relative
  //  lcp (int): bra orbital angular momentum, cm
  //  lc  (int): ket orbital angular momentum, cm
  //  Lp  (int): bra total orbital angular momentum
  //  L   (int): ket total orbital angular momentum
  //  Sp  (int): bra spin
  //  S   (int): ket spin
  //  Jp  (int): bra angular momentum
  //  J   (int): ket angular momentum
  // Returns:
  //  reduced matrix element (double)
  inline double RelativeCMLsumRME(const int& lrp, const int& lr,
                                  const int& lcp, const int& lc,
                                  const int& Lp, const int& L,
                                  const int& Sp, const int& S,
                                  const int& Jp, const int& J)
  {
    if (lrp != lr && lcp != lc)
      return 0;
    auto hat_product = HatProduct(Lp, J, lr, lc, L);
    auto parity = ParitySign(J + Lp + S + lc + lr);
    auto LS6j = am::Wigner6J(L, Lp, 1, Jp, J, S);
    auto cm_term = (ParitySign(Lp) * am::Wigner6J(lc, lc, 1, L, Lp, lr)
                    * (AngularMomentumJRME(lc, lc) / Hat(lr)));
    auto rel_term = (ParitySign(L) * am::Wigner6J(lr, lr, 1, Lp, L, lc)
                     * (AngularMomentumJRME(lr, lr) / Hat(lc)));
    auto result = (hat_product * parity * LS6j * (cm_term + rel_term));
    return result;
  }

  // Calculate matrix element of the gradient operator in radial HO basis.
  // The HO length parameter is taken to be 1.
  //
  // Arguments:
  //  np (int): bra radial quantum number
  //  n  (int): ket radial quantum number
  //  lp (int): bra orbital angular momentum
  //  l  (int): ket orbital angular momentum
  // Returns:
  //  matrix element (double)
  inline double GradientME(const int& np, const int& n,
                           const int& lp, const int& l)
  {
    double result = 0;
    if (lp == l + 1)
      {
        result = (-std::sqrt(n + l + 1.5) * (np == n)
                  + 3 * std::sqrt(n) * (np == n - 1));
        result *= std::sqrt((l + 1.0) / (2 * l + 3.0));
      }
    if (l == lp + 1)
      {
        result = (-std::sqrt(np + lp + 1.5) * (n == np)
                  + 3 * std::sqrt(np) * (n == np - 1));
        result *= std::sqrt((lp + 1.0) / (2 * lp + 1));
      }
    return result;
  }

  // Calculate matrix element of r (or R) in radial HO basis.
  // The HO length parameter is taken to be 1.
  //
  // Arguments:
  //  np (int): bra radial quantum number
  //  n  (int): ket radial quantum number
  //  lp (int): bra orbital angular momentum
  //  l  (int): ket orbital angular momentum
  // Returns:
  //  matrix element (double)
  inline double RadiusME(const int& np, const int& n,
                         const int& lp, const int& l)
  {
    double result = 0;
    if (lp == l + 1)
      {
        result = (std::sqrt(n + l + 1.5) * (np == n)
                  - std::sqrt(n) * (np + 1 == n));
      }
    if (l == lp + 1)
      {
        result = RadiusME(n, np, l, lp);
      }
    if (lp == l)
      {
        quadrature::gsl_params_2n p{np, lp, n, l, false, 1, 1};
        auto norm_product = (ho::CoordinateSpaceNorm(n, l, 1)
                            * ho::CoordinateSpaceNorm(np, lp, 1));
        result = norm_product * quadrature::IntegralMPiR(p);
      }
    result *= SphericalHarmonicCRME(lp, l, 1);
    return result;
  }
}

#endif
