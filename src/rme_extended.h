/*
 * This adds to the reduced matrix elements in rme.h of the am library. The
 * reduced matrix elements are defined using the Rose convention,
 * <j' m'|T_{LM}|j m> = (j m L M|j' m') <j'||T_L||j>.
 */

#ifndef RME_EXTENDED_H
#define RME_EXTENDED_H

#include "am/rme.h"
#include "constants.h"

namespace am
{
  // Calculates reduced matrix element of the dot product of the Pauli
  // matrices acting on nucleons 1 and 2.
  //
  // Arguments:
  //   sp(int): Final (iso)spin
  //   s(int) : Initial (iso)spin
  // Returns:
  //   reduced matrix element (double), Rose convention
  inline double PauliDotProductRME(const int& sp, const int& s)
  {
    double result = (-3 * (s == 0) + (s == 1)) * (sp == s);
    return result;
  }

  // Calculates reduced matrix element of the total spin operator \vec{S} of a
  // two nucleon system in the LS-coupling scheme.
  //
  // Arguments:
  //   lp(int): Final orbital angular momentum
  //   sp(int): Final spin
  //   jp(int): Final total angular momentum
  //   l(int) : Initial orbital angular momentum
  //   s(int) : Initial spin
  //   j(int) : Initial total angular momentum
  // Returns:
  //   reduced matrix element in coupled LS basis (double), Rose convention
  inline double LSCoupledTotalSpinRME(const int& lp, const int& sp, const int& jp,
                                      const int& l, const int& s, const int& j)
  {
    bool kronecker = (lp == l && sp == 1 && s == 1);
    if (!kronecker)
      return 0;
    auto result = std::sqrt(6) * std::pow(-1, jp + l) * Hat(j);
    result *= am::Wigner6J(1, jp, l, j, 1, 1);
    return result;
  }

  // Calculates reduced matrix element of [\vec{S} ⊗ Y_0(\hat{r})]_1.
  //
  // Arguments:
  //   lp(int): Final orbital angular momentum
  //   sp(int): Final spin
  //   jp(int): Final total angular momentum
  //   l(int) : Initial orbital angular momentum
  //   s(int) : Initial spin
  //   j(int) : Initial total angular momentum
  // Returns:
  //   reduced matrix element in coupled LS basis (double), Rose convention
  inline double LSCoupledTotalSpinY0Rank1RME(const int& lp, const int& sp, const int& jp,
                                             const int& l, const int& s, const int& j)
  {
    bool kronecker = (lp == l && sp == 1 && s == 1);
    if (!kronecker)
      return 0;
    auto result = Hat(j) * std::pow(-1, l + jp) * am::Wigner6J(1, j, l, jp, 1, 1);
    result *= std::sqrt(3.0 / (2.0 * constants::pi));
    return result;
  }

  // Calculates reduced matrix element of [\vec{S} ⊗ Y_2(\hat{r})]_1.
  //
  // Arguments:
  //   lp(int): Final orbital angular momentum
  //   sp(int): Final spin
  //   jp(int): Final total angular momentum
  //   l(int) : Initial orbital angular momentum
  //   s(int) : Initial spin
  //   j(int) : Initial total angular momentum
  // Returns:
  //   reduced matrix element in coupled LS basis (double), Rose convention
  inline double LSCoupledTotalSpinY2Rank1RME(const int& lp, const int& sp, const int& jp,
                                             const int& l, const int& s, const int& j)
  {
    auto result = Hat(j) * Hat(l) * Hat(lp) * std::pow(-1, l);
    result *= am::Wigner9J(l, j, 1, lp, jp, 1, 2, 1, 1);
    result *= am::Wigner3J(l, 2, lp, 0, 0, 0);
    result *= 3 * std::sqrt(5 / (2 * constants::pi));
    return result;
  }
}

#endif
