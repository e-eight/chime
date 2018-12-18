#ifndef CHIRAL_H
#define CHIRAL_H

#include "basis/lsjt_scheme.h"

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
        enum class Name
        {
            identity,
            charge_radius,
            gamow_teller,
        };
        // Chiral order of operators
        enum class Order { lo, nlo, n2lo, n3lo, n4lo, full };
        
        // Constructors
        Operator() {}

        // Accessors
        double get_rme() const { return rme; }
        void set_rme()
        {
            switch(name)
            {
            case chiral::Operator::Name::charge_radius:
                rme = chiral::operator_rc_sq(order, bra, ket, osc_b);
                break;
            case chiral::Operator::Name::gamow_teller:
                rme = chiral::operator_gt(order, bra, ket, osc_b);
                break;
            default:
                rme = 1.0;
                break;
            }
        }
        
        // Destructor
        ~Operator() {}

        // Public members
        Operator::Name name;
        Operator::Order order;
        basis::RelativeStateLSJT bra;
        basis::RelativeStateLSJT ket;
        double osc_b; // Oscillator constant b

    private:
        double rme; // Reduced matrix element
    };
}

#endif
