#include "chiral.h"

namespace chiral
{
  double Operator::ReducedMatrixElement(const Order order,
                                        const basis::RelativeStateLSJT& bra,
                                        const basis::RelativeStateLSJT& ket,
                                        const double& osc_b)
  {
    switch(order)
      {
      case Order::lo:
        return LOMatrixElement(bra, ket, osc_b);
      case Order::nlo:
        return NLOMatrixElement(bra, ket, osc_b);
      case Order::n2lo:
        return N2LOMatrixElement(bra, ket, osc_b);
      case Order::n3lo:
        return N3LOMatrixElement(bra, ket, osc_b);
      case Order::n4lo:
        return N4LOMatrixElement(bra, ket, osc_b);
      }
  }
}
