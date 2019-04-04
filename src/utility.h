#ifndef UTILITY_H
#define UTILITY_H

#include <memory>

namespace util
{
  inline int KroneckerDelta(int a, int b) { return a == b; };

  inline double PauliDotProduct(int s)
  {
    // Returns the matrix element of \sigma_1 \dot \sigma_2.
    if (s == 0)
      return -3.0;
    if (s == 1)
      return 1.0;
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
