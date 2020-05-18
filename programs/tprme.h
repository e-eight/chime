/*******************************************************************************
 tprme.h

 Defines functions for calculating the reduced matrix elements of some LS
 coupled tensor products, in the Rose convention.

 Language: C++11
 Soham Pal
 Iowa State University
*******************************************************************************/

#ifndef TPRME_H_
#define TPRME_H_

#include <cmath>

#include "am/rme.h"
#include "basis/lsjt_scheme.h"

namespace chime {
namespace tp {

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

// Calculate reduced matrix element of tensor products of Pauli matrices.
// Possible ranks of the tensor products are 0, 1, 2.
//
// Arguments:
//   sp (int): bra spin of two nucleon system
//   s (int): ket spin of two nucleon system
//   rank (int): rank of tensor product of Pauli matrices
// Returns:
//   reduced matrix element
inline double SpinTensorProductRME(const int& sp, const int& s, const int& rank)
{
  assert((rank >= 0) && (rank <= 2));
  assert((sp == 1) || (sp == 0) || (s == 1) || (s == 0));

  if (rank == 0) {
    return ParitySign(s) * Hat(1 - s) * (1 - std::abs(sp - s)) / Hat(s);
  }
  else if (rank == 1) {
    return std::sqrt(2) * Hat(s) * std::abs(sp - s);
  }
  else {
    return std::sqrt(20. / 3.) * (1 - std::abs(sp - s))
           * (1 - std::abs(sp - 1));
  }
}

// Calculate reduced matrix element of rank `c` tensor product of C_a and
// Σ_b, where C_a is a rank `a` spherical harmonic with Racah normalization,
// and Σ_b is a rank `b` tensor product of Pauli matrices, in the LS coupled
// scheme.
//
// Arguments:
//   bra (basis::RelativeSubspaceLSJT): bra subspace containing am labels
//   ket (basis::RelativeSubspaceLSJT): ket subspace containing am labels
//   a (int): rank of spherical harmonic
//   b (int): rank of Pauli matrix tensor product
//   c (int): rank of operator
// Returns:
//   reduced matrix element
inline double CSpinTensorProductRME(const basis::RelativeSubspaceLSJT& bra,
                                    const basis::RelativeSubspaceLSJT& ket,
                                    const int& a, const int& b, const int& c)
{
  // Get angular momentum labels.
  int bra_L = bra.L();
  int bra_S = bra.S();
  int bra_J = bra.J();
  int ket_L = ket.L();
  int ket_S = ket.S();
  int ket_J = ket.J();

  if (am::AllowedTriangle(bra_L, a, ket_L)) {
    double result = HatProduct(bra_L, bra_S, ket_J, c);
    result *= am::Wigner9J(ket_L, ket_S, ket_J, a, b, c, bra_L, bra_S, bra_J);
    result *= am::SphericalHarmonicCRME(bra_L, ket_L, a);
    result *= SpinTensorProductRME(bra_S, ket_S, b);
    return result;
  }
  return 0;
}

// Calculate reduced matrix element of rank `e` tensor product of
// [C_a C_b]_c and Σ_d in the llLS coupled scheme. Here [C_a C_b]_c is a
// rank `c` tensor product of C_a and C_b; C_a and C_b are spherical
// harmonics with Racah normalization, of ranks `a` and `b`, respectively.
// C_a and C_b are associated with different coordinates. For the use cases
// in `chime`, these are the relative and cm coordinates, respectively.
// Σ_d is a rank `d` tensor product of Pauli matrices.
//
// Arguments:
//   bra (basis::RelativeCMStateLSJT): bra state containing am labels
//   ket (basis::RelativeCMStateLSJT): ket state containing am labels
//   a (int): rank of relative spherical harmonic
//   b (int): rank of cm spherical harmonic
//   c (int): rank of spherical harmonic tensor product
//   d (int): rank of Pauli matrix tensor product
//   e (int): rank of operator
// Returns:
//   reduced matrix element
inline double CCSpinTensorProductRME(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const int& a, const int& b, const int& c,
                                     const int& d, const int& e)
{
  // Get angular momentum labels.
  int bra_lr = bra.lr();
  int bra_lc = bra.lc();
  int bra_L = bra.L();
  int bra_S = bra.S();
  int bra_J = bra.J();
  int ket_lr = ket.lr();
  int ket_lc = ket.lc();
  int ket_L = ket.L();
  int ket_S = ket.S();
  int ket_J = ket.J();

  if (am::AllowedTriangle(bra_lr, a, ket_lr)
      && am::AllowedTriangle(bra_lc, b, ket_lc)) {
    double result =
        HatProduct(bra_L, bra_S, ket_J, e, bra_lr, bra_lc, ket_L, c);
    result *= am::Wigner9J(ket_L, ket_S, ket_J, c, d, e, bra_L, bra_S, bra_J);
    result *=
        am::Wigner9J(ket_lr, ket_lc, ket_L, a, b, c, bra_lr, bra_lc, bra_L);
    result *= am::SphericalHarmonicCRME(bra_lr, ket_lr, a);
    result *= am::SphericalHarmonicCRME(bra_lc, ket_lc, b);
    result *= SpinTensorProductRME(bra_S, ket_S, d);
    return result;
  }
  return 0;
}
}  // namespace tp
}  // namespace chime

#endif
