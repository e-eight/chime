#ifndef UTILITY_H
#define UTILITY_H

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

#endif
