#ifndef CHIRAL_H
#define CHIRAL_H

#include "basis/lsjt_scheme.h"

/* 

@file chiral.h

Author(s): Soham Pal
Iowa State University

*/

namespace chiral
{
    // Operator Names
    enum class Name
    {
        identity,
        charge_radius,
        gamow_teller
    };

    // Chiral Orders
    enum class Order
    { lo, nlo, n2lo, n3lo, n4lo, full };

    class ChiralOperator
    {
    public:
        // Constructors
        ChiralOperator();
        ChiralOperator(Name name);
        ChiralOperator(Name name, Order order);
        ChiralOperator(Name name, int J0, int T0);
        ChiralOperator(Name name, Order order, int J0, int T0);
        
        // Destructor
        virtual ~ChiralOperator() = 0;

        // Methods
        virtual void calculate_rme(const basis::RelativeStateLSJT& bra,
                                   const basis::RelativeStateLSJT& ket,
                                   const double& osc_b,
                                   double& rme) = 0;
        // Data Members 
        Name name;
        Order order;
        int J0; // Tensor rank
        int T0; // Isotensor rank
    };
}

#endif
