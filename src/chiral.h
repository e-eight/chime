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
    struct ChiralOperator
    {
        // Operator Names
        enum struct Name
        {
            identity,
            charge_radius
            // gamow_teller
        };

        // Chiral Orders
        enum struct Order
        { lo, nlo, n2lo, n3lo, n4lo, full };

        ChiralOperator::Name name = ChiralOperator::Name::identity;
        ChiralOperator::Order order = ChiralOperator::Order::full;
    };

    void calculate_rme(const ChiralOperator& op,
                       const basis::RelativeStateLSJT& bra,
                       const basis::RelativeStateLSJT& ket,
                       const double& osc_b,
                       double& rme);
}

#endif
