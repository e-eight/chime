#ifndef CHIRAL_H
#define CHIRAL_H

#include "basis/lsjt_scheme.h"
#include "enum.h"
#include "factory.h"

namespace chiral
{
  // Chiral Orders
  BETTER_ENUM (Order, int, LO=1, NLO, N2LO, N3LO, N4LO);

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
                                   const double& osc_b,
                                   const double& regulator) = 0;
    virtual double LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                   const basis::RelativeCMStateLSJT& ket,
                                   const double& osc_b,
                                   const double& regulator) = 0;

    virtual double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                    const basis::RelativeStateLSJT& ket,
                                    const double& osc_b,
                                    const double& regulator) = 0;
    virtual double NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                    const basis::RelativeCMStateLSJT& ket,
                                    const double& osc_b,
                                    const double& regulator) = 0;

    virtual double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator) = 0;
    virtual double N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator) = 0;

    virtual double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator) = 0;
    virtual double N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator) = 0;

    virtual double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator) = 0;
    virtual double N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const double& osc_b,
                                     const double& regulator) = 0;

    /* double RelativeRME(const Order& order, */
    /*                    const basis::RelativeStateLSJT& bra, */
    /*                    const basis::RelativeStateLSJT& ket, */
    /*                    const double& osc_b, */
    /*                    const double& regulator); */

    /* double RelativeCMRME(const Order& order, */
    /*                      const basis::RelativeCMStateLSJT& bra, */
    /*                      const basis::RelativeCMStateLSJT& ket, */
    /*                      const double& osc_b, */
    /*                      const double& regulator); */
  };

  template <class StateType, class OscillatorType>
    double ReducedMatrixElement(const Operator& op,
                                const Order& ord,
                                const StateType& bra,
                                const StateType& ket,
                                const OscillatorType& b)
  {
    switch(ord)
      {
      case Order::LO:
        return op->LOMatrixElement(bra, ket, b);
      case Order::NLO:
        return op->NLOMatrixElement(bra, ket, b);
      case Order::N2LO:
        return op->N2LOMatrixElement(bra, ket, b);
      case Order::N3LO:
        return op->N3LOMatrixElement(bra, ket, b);
      case Order::N4LO:
        return op->N4LOMatrixElement(bra, ket, b);
      }
  }


}

#endif
