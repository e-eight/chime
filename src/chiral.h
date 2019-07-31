#ifndef CHIRAL_H
#define CHIRAL_H

#include "basis/lsjt_scheme.h"
#include "enum.h"
#include "utility.h"
#include "factory.h"

namespace chiral
{
  // Chiral orders
  BETTER_ENUM(Order, int, LO=1, NLO, N2LO, N3LO, N4LO);

  // Chiral operator interface
  struct Operator : factory::Factory<Operator>
  {
    Operator(Key) {}

    virtual ~Operator() = default;

    virtual int G0() = 0; // returns parity of operator

    virtual int J0() = 0; // returns tensor rank of operator

    virtual int T0() = 0; // returns isotensor rank of operator

    virtual double LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                   const basis::RelativeStateLSJT& ket,
                                   const util::OscillatorParameter& b,
                                   const bool& regularize,
                                   const double& regulator) = 0;
    virtual double LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                   const basis::RelativeCMStateLSJT& ket,
                                   const util::OscillatorParameter& b,
                                   const bool& regularize,
                                   const double& regulator) = 0;

    virtual double NLOMatrixElement(const basis::RelativeStateLSJT& bra,
                                    const basis::RelativeStateLSJT& ket,
                                    const util::OscillatorParameter& b,
                                    const bool& regularize,
                                    const double& regulator) = 0;
    virtual double NLOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                    const basis::RelativeCMStateLSJT& ket,
                                    const util::OscillatorParameter& b,
                                    const bool& regularize,
                                    const double& regulator) = 0;

    virtual double N2LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator) = 0;
    virtual double N2LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator) = 0;

    virtual double N3LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator) = 0;
    virtual double N3LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator) = 0;

    virtual double N4LOMatrixElement(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator) = 0;
    virtual double N4LOMatrixElement(const basis::RelativeCMStateLSJT& bra,
                                     const basis::RelativeCMStateLSJT& ket,
                                     const util::OscillatorParameter& b,
                                     const bool& regularize,
                                     const double& regulator) = 0;

  template <class StateType>
  double ReducedMatrixElement(const Order& ord,
                              const StateType& bra,
                              const StateType& ket,
                              const util::OscillatorParameter& b,
                              const double& regulator,
                              const bool& regularize)
    {
    switch(ord)
      {
      case Order::LO:
        return LOMatrixElement(bra, ket, b, regulator, regularize);
      case Order::NLO:
        return NLOMatrixElement(bra, ket, b, regulator, regularize);
      case Order::N2LO:
        return N2LOMatrixElement(bra, ket, b, regulator, regularize);
      case Order::N3LO:
        return N3LOMatrixElement(bra, ket, b, regulator, regularize);
      case Order::N4LO:
        return N4LOMatrixElement(bra, ket, b, regulator, regularize);
      }
  }


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

}

#endif
