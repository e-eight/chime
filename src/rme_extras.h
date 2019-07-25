/*
 * Defines reduced matrix elements for two-body Pauli matrices,
 * [\hat{r} ⊗ \hat{r}]_λ, and their tensor products.
 */

#ifndef RME_EXTRAS_H
#define RME_EXTRAS_H

#include <cmath>
#include "am/wigner_gsl.h"
#include "am/racah_reduction.h"

namespace  am
{
  enum class Nucleon : int { one = 1, two = 2 };

  inline double ssCoupledPauliRME(const int& sp, const int& s, const Nucleon& i)
  {
    // Calculate reduced matrix element of Pauli matrix in a two-body coupled
    // (iso)spin basis.
    //
    // Arguments:
    //  sp (int): bra spin
    //  s (int): ket spin
    //  i (Nucleon): nucleon index;
    // Returns:
    //  reduced matrix element (double), Rose convention

    assert(AllowedTriangle(HalfInt(1, 2), HalfInt(1, 2), sp));
    assert(AllowedTriangle(HalfInt(1, 2), HalfInt(1, 2), s));

    double result = 0;
    if (sp == s)
      result = std::sqrt(s * (s + 1.));
    else if (i == Nucleon::one)
      result = ParitySign(s) * std::sqrt(3);
    else if (i == Nucleon::two)
      result = ParitySign(sp) * std::sqrt(3);

    return result;
  }

  inline double LSCoupledPauliRME(const int& lp,
                                  const int& sp,
                                  const int& jp,
                                  const int& l,
                                  const int& s,
                                  const int& j,
                                  const Nucleon& i)
  {
    // Calculate reduced matrix element of Pauli matrix in a LS-coupled spin
    // basis.
    //
    // Arguments:
    //  lp (int): bra orbital angular momentum
    //  sp (int): bra spin
    //  jp (int): bra total angular momentum
    //  l (int): ket orbital angular momentum
    //  s (int): ket spin
    //  j (int): ket total angular momentum
    //  k (int): rank of [\vec{σ}_1 ⊗ \vec{σ}_2]
    // Returns:
    //  reduced matrix element (double), Rose convention
    if (lp != l)
      return 0;
    double hat_product = Hat(j) * Hat(sp);
    double wigner_6j = ParitySign(s+l+jp+1) * Wigner6J(s, j, l, jp, sp, 1);
    double rme_product = ssCoupledPauliRME(sp, s, i);
    double result = hat_product * wigner_6j * rme_product;

    return result;
  }

  inline double ssCoupledPauliProductRME(const int& sp, const int& s, const int& k)
  {
    // Calculate reduced matrix element of [\vec{σ}_1 ⊗ \vec{σ}_2]_k in a
    // two-body coupled (iso)spin basis.
    //
    // Arguments:
    //  sp (int): bra spin
    //  s (int): ket spin
    //  k (int): rank of the tensor product;
    // Returns:
    //  reduced matrix element (double), Rose convention

    assert(AllowedTriangle(1, 1, k));

    double hat_product = 2 * Hat(s) * Hat(k);
    HalfInt half = HalfInt(1, 2);
    double wigner_9j = Wigner9J(half, half, s, 1, 1, k, half, half, sp);
    double result = 3 * hat_product * wigner_9j;

    return result;
  }

  inline double LSCoupledPauliProductRME(const int& lp,
                                         const int& sp,
                                         const int& jp,
                                         const int& l,
                                         const int& s,
                                         const int& j,
                                         const int& k)
  {
    // Calculate reduced matrix element of [\vec{σ}_1 ⊗ \vec{σ}_2]_k in a
    // LS-coupled basis.
    //
    // Arguments:
    //  lp (int): bra orbital angular momentum
    //  sp (int): bra spin
    //  jp (int): bra total angular momentum
    //  l (int): ket orbital angular momentum
    //  s (int): ket spin
    //  j (int): ket total angular momentum
    //  k (int): rank of [\vec{σ}_1 ⊗ \vec{σ}_2]
    // Returns:
    //  reduced matrix element (double), Rose convention
    if (lp != l)
      return 0;
    double hat_product = Hat(j) * Hat(sp);
    double wigner_6j = ParitySign(s+l+jp+1) * Wigner6J(s, j, l, jp, sp, k);
    double rme_product = ssCoupledPauliProductRME(sp, s, k);
    double result = hat_product * wigner_6j * rme_product;

    return result;
  }

  inline double RHatProductRME(const int& lp, const int& l, const int& k)
  {
    // Calculate reduced matrix element of [\hat{r} ⊗ \hat{r}]_k.
    //
    // Arguments:
    //  lp (int): bra orbital angular momentum
    //  l (int): ket orbital angular momentum
    //  k (int): rank of [\hat{r} ⊗ \hat{r}]
    // Returns:
    //  reduced matrix element (double), Rose convention

    assert(AllowedTriangle(1, 1, k));

    // parity constraint
    if ((lp + l + k) % 2 != 0)
      return 0;

    double hat_product = ParitySign(lp + k) * Hat(l) * Hat(k);
    double wigner_3j_product = (Wigner3J(1, k, 1, 0, 0, 0)
                                * Wigner3J(lp, k, l, 0, 0, 0));
    double result = hat_product * wigner_3j_product;
    return result;
  }

  inline double LSCoupledRHatProductPauliRME(const int& lp,
                                             const int& sp,
                                             const int& jp,
                                             const int& l,
                                             const int& s,
                                             const int& j,
                                             const int& k1,
                                             const int& k2,
                                             const Nucleon& i)
  {
    // Calculate reduced matrix element of [[\hat{r} ⊗ \hat{r}]_k1 ⊗ \vec{σ}_i]_k2
    // in a LS-coupled two body basis.
    //
    // Arguments:
    //  lp (int): bra orbital angular momentum
    //  sp (int): bra spin
    //  jp (int): bra total angular momentum
    //  l (int): ket orbital angular momentum
    //  s (int): ket spin
    //  j (int): ket total angular momentum
    //  k1 (int): rank of [\hat{r} ⊗ \hat{r}]
    //  k2 (int): rank of [[\hat{r} ⊗ \hat{r}]_k1 ⊗ \vec{σ}_i]
    //  i (Nucleon): nucleon index
    // Returns:
    //  reduced matrix element (double), Rose convention

    assert(AllowedTriangle(k1, 1, k2));

    double hat_product = Hat(j) * Hat(k2) * Hat(lp) * Hat(sp);
    double wigner_9j = Wigner9J(l, s, j, k1, 1, k2, lp, sp, jp);
    double rme_product = RHatProductRME(lp, l, k1) * ssCoupledPauliRME(sp, s, i);
    double result = hat_product * wigner_9j * rme_product;

    return result;
  }

  inline double LSCoupledRHatProductPauliProductRME(const int& lp,
                                                    const int& sp,
                                                    const int& jp,
                                                    const int& l,
                                                    const int& s,
                                                    const int& j,
                                                    const int& k1,
                                                    const int& k2,
                                                    const int& k3)
  {
    // Calculate reduced matrix element of
    // [[\hat{r} ⊗ \hat{r}]_k1 ⊗ [\vec{σ}_1 ⊗ \vec{σ}_2]_k2]_k3
    // in a LS-coupled two body basis.
    //
    // Arguments:
    //  lp (int): bra orbital angular momentum
    //  sp (int): bra spin
    //  jp (int): bra total angular momentum
    //  l (int): ket orbital angular momentum
    //  s (int): ket spin
    //  j (int): ket total angular momentum
    //  k1 (int): rank of [\hat{r} ⊗ \hat{r}];
    //  k2 (int): rank of [\vec{σ}_1 ⊗ \vec{σ}_2]
    //  k3 (int): rank of [[\hat{r} ⊗ \hat{r}]_k1 ⊗ [\vec{σ}_1 ⊗ \vec{σ}_2]_k2]
    // Returns:
    //  reduced matrix element (double), Rose convention

    assert(AllowedTriangle(k1, k2, k3));

    double hat_product = Hat(j) * Hat(k3) * Hat(lp) * Hat(sp);
    double wigner_9j = Wigner9J(l, s, j, k1, k2, k3, lp, sp, jp);
    double rme_product = (RHatProductRME(lp, l, k1)
                          * ssCoupledPauliProductRME(sp, s, k2));
    double result = hat_product * wigner_9j * rme_product;

    return result;
  }
}

#endif
