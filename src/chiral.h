#ifndef CHIRAL_H
#define CHIRAL_H

#include "lib/basis/lsjt_scheme.h"

/* 

@file chiral.h

Defines the base chiral operator class for the computation of the matrix
elements.

The class members are:
- operator name (all implemented operators are listed in enum class Names)
- chiral order (all implemented operators are listed in enum class orders)
- bra and ket defining the final and initial states (depends on ND basis)
- oscillator constant
- reduced matrix element

Any new chiral operator must first be added to the enum class Names, and
similarly any new chiral order must first be added to the enum class Orders.

Author(s): Soham Pal
Iowa State University

*/

namespace chiral
{
    class Operator
    {
    public:
        // Operator names
        enum class Names
        {
            identity,
            charge_radius,
            gamow_teller,
        };
        // Chiral order of operators
        enum class Orders { lo, nlo, n2lo, n3lo, n4lo, full };
        
        // Constructors
        Operator() {}
        Operator(const Operator::Names name): _name{name} {}

        // Accessors
        Operator::Names get_name() const { return _name; }
        void set_order(const Operator::Orders order) { _order = order; }
        void set_bra(const basis::RelativeStateLSJT& bra) { _bra = bra; }
        void set_ket(const basis::RelativeStateLSJT& ket) { _ket = ket; }
        void set_b(const double b) { osc_constant = b; }
        virtual void set_rme() { rme = 1.0; }
        double get_rme() const { return rme; }

        // Destructor
        virtual ~Operator() {}

    protected:
        Operator::Names _name;
        Operator::Orders _order;
        basis::RelativeStateLSJT _bra;
        basis::RelativeStateLSJT _ket;
        double osc_constant; // Oscillator constant b
        double rme; // Reduced matrix element
    };
}

#endif
