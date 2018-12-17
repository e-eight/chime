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
        double get_rme() const { return rme; }
        virtual void set_rme() { rme = 1.0; }

        // Destructor
        virtual ~Operator() {}

        // Public members
        Operator::Orders _order;
        basis::RelativeStateLSJT _bra;
        basis::RelativeStateLSJT _ket;
        double osc_constant; // Oscillator constant b

    private:
        Operator::Names _name;
        double rme; // Reduced matrix element
    };
}

#endif
