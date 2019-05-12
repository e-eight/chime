/*
 * This adds to the reduced matrix elements in rme.h of the am library. The
 * reduced matrix elements are defined using the Rose convention,
 * <j' m'|T_{LM}|j m> = (j m L M|j' m') <j'||T_L||j>.
 */

#ifndef RME_EXTENDED_H
#define RME_EXTENDED_H

#include "basis/am/rme.h"

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

  // Calculates reduced matrix element of the total spin operator of a two
  // nucleon system in the LS-coupling scheme.
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
    bool kronecker = lp == l && sp == 1 && s == 1;
    if (!kronecker)
      return 0;
    auto result = std::sqrt(6) * std::pow(-1, jp + l) * Hat(j);
    result *= am::Wigner6J(1, jp, l, j, 1, 1);
    return result;
  }

  // Calculates reduced matrix element of \vec{S}.\hat{r}\hat{r}, where \hat{r}
  // is the unit vector in the radial direction, and \vec{S} is the total spin of
  // the two nucleon system, in the LS-coupling scheme.
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
  inline double LSCoupledSDotRhatRhatRME(const int& lp, const int& sp, const int& jp,
                                         const int& l, const int& s, const int& j)
  {
    bool kronecker = sp == 1 && s == 1 && am::AllowedTriangle(j, jp, 1);
    if (!kronecker)
      return 0;
    double result = 0;
    for (int k = std::abs(l - 1); k <= l + 1; ++k)
      {
        auto product6j = (am::Wigner6J(lp, 1, jp, 1, k, 1)
                          * am::Wigner6J(1, k, jp, 1, j, l));
        auto product3j = (am::Wigner3J(lp, 1, k, 0, 0, 0)
                          * am::Wigner3J(k, 1, l, 0, 0, 0));
        result += Hat(k) * Hat(k) * product6j * product3j;
      }
    result *= Hat(j) * Hat(lp) * Hat(l);
    result *= std::sqrt(6) * std::pow(-1, jp + j);
    return result;
  }
}

#endif
