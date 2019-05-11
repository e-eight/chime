#ifndef UTILITY_H
#define UTILITY_H

#include <memory>
#include <cmath>
#include "basis/am/halfint.h"
#include "basis/am/wigner_gsl.h"

namespace util
{
  inline int KroneckerDelta(int a, int b) { return a == b; };

  inline double PauliDotProduct(int si, int sf)
  {
    // Returns the reduced matrix element of \sigma_1 \cdot \sigma_2.

    return (-3 * (si == 0) + (si == 1)) * (si == sf);
  }

  inline double TotalSpin(int li, int si, int ji, int lf, int sf, int jf)
  {
    // Returns the reduced matrix element of S = \sigma_1 + \sigma_2 in the
    // |lsj> basis.
    if (si == 1 && sf == 1)
      return (std::sqrt(6) * std::pow(-1, jf + li) * Hat(ji)
              * am::Wigner6J(1, jf, li, ji, 1, 1) * (li == lf));
    return 0;
  }

  // make_unique, because C++11 does not have one.
  // https://herbsutter.com/gotw/_102/ for details.
  template <class T, class ...Args>
  std::unique_ptr<T> make_unique(Args&&... args)
  {
    return std::unique_ptr<T>( new T( std::forward<Args>(args)... ));
  }
}
#endif
