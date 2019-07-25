#include "chiral.h"

namespace chiral
{
  double Operator::RelativeRME(const Order& order,
                               const basis::RelativeStateLSJT& bra,
                               const basis::RelativeStateLSJT& ket,
                               const double& osc_b,
                               const double& regulator)
  {
    switch(order)
      {
      case Order::LO:
        return LOMatrixElement(bra, ket, osc_b, regulator);
      case Order::NLO:
        return NLOMatrixElement(bra, ket, osc_b, regulator);
      case Order::N2LO:
        return N2LOMatrixElement(bra, ket, osc_b, regulator);
      case Order::N3LO:
        return N3LOMatrixElement(bra, ket, osc_b, regulator);
      case Order::N4LO:
        return N4LOMatrixElement(bra, ket, osc_b, regulator);b
      }
  }

  double Operator::RelativeCMRME(const Order& order,
                                 const basis::RelativeCMStateLSJT& bra,
                                 const basis::RelativeCMStateLSJT& ket,
                                 const double& osc_b,
                                 const double& regulator)
  {
    switch(order)
      {
      case Order::LO:
        return LOMatrixElement(bra, ket, osc_b, regulator);
      case Order::NLO:
        return NLOMatrixElement(bra, ket, osc_b, regulator);
      case Order::N2LO:
        return N2LOMatrixElement(bra, ket, osc_b, regulator);
      case Order::N3LO:
        return N3LOMatrixElement(bra, ket, osc_b, regulator);
      case Order::N4LO:
        return N4LOMatrixElement(bra, ket, osc_b, regulator);
      }
  }
}
