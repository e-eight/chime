#ifndef CHIRAL_H
#define CHIRAL_H

#include "basis/lsjt_scheme.h"
#include "factory.h"

namespace chiral
{
  // Chiral Orders
  enum struct Order { lo, nlo, n2lo, n3lo, n4lo };

  // Chiral Operator
  struct Operator : factory::Factory<Operator>
  {
    Operator(Key) {}

    virtual ~Operator() = default;

    virtual int G0() = 0; // returns parity of operator

    virtual int J0() = 0; // returns tensor rank of operator

    virtual int T0() = 0; // returns isotensor rank of operator

    virtual double LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                   const basis::RelativeStateLSJT& ket,
                                   const double& osc_b) = 0;

    virtual double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                    const basis::RelativeStateLSJT& ket,
                                    const double& osc_b) = 0;

    virtual double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b) = 0;

    virtual double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b) = 0;

    virtual double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b) = 0;

    double ReducedMatrixElement(const Order order,
                                const basis::RelativeStateLSJT& bra,
                                const basis::RelativeStateLSJT& ket,
                                const double& osc_b);
  };
}

#endif
