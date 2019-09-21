#ifndef CHIRAL_H
#define CHIRAL_H

#include <vector>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <stdexcept>
#include "basis/lsjt_scheme.h"
#include "utility.h"
#include "factory.h"

namespace chiral
{
  // Chiral orders
  enum struct Order {lo, nlo, n2lo, n3lo, n4lo};

  const std::vector<Order> v_order
    {Order::lo, Order::nlo, Order::n2lo, Order::n3lo, Order::n4lo};

  static std::unordered_map<std::string, Order> m_order
    {
     {"lo", Order::lo},
     {"nlo", Order::nlo},
     {"n2lo", Order::n2lo},
     {"n3lo", Order::n3lo},
     {"n4lo", Order::n4lo}
    };

  static std::unordered_map<Order, std::string> reverse_m_order
    {
     {Order::lo, "lo"},
     {Order::nlo, "nlo"},
     {Order::n2lo, "n2lo"},
     {Order::n3lo, "n3lo"},
     {Order::n4lo, "n4lo"}
    };

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
                                const bool& regularize,
                                const double& regulator)
    {
      auto found = std::find(v_order.begin(), v_order.end(), ord);
      if (found == v_order.end())
        throw std::invalid_argument("Chiral order must be in [lo, n4lo]");
      switch(ord)
        {
        case Order::lo:
          return LOMatrixElement(bra, ket, b, regularize, regulator);
        case Order::nlo:
          return NLOMatrixElement(bra, ket, b, regularize, regulator);
        case Order::n2lo:
          return N2LOMatrixElement(bra, ket, b, regularize, regulator);
        case Order::n3lo:
          return N3LOMatrixElement(bra, ket, b, regularize, regulator);
        case Order::n4lo:
          return N4LOMatrixElement(bra, ket, b, regularize, regulator);
        }
    }
  };

}

#endif
