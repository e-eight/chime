#ifndef CHIRAL_H
#define CHIRAL_H

#include <string>
#include <memory>
#include "basis/lsjt_scheme.h"

/* 

@file chiral.h

Author(s): Soham Pal
Iowa State University

*/

namespace chiral
{
    // Chiral Orders
    enum class Order
    { lo, nlo, n2lo, n3lo, n4lo, full };

    class ChiralOperator
    {
    public:
        // Constructors
        ChiralOperator();
        ChiralOperator(int J0, int T0);
        ChiralOperator(int G0, int J0, int T0);
        
        // Destructor
        virtual ~ChiralOperator() = 0;

        // Methods
        using operator_ptr = std::unique_ptr<ChiralOperator>;
        static operator_ptr create_operator(std::string name);
        virtual double calculate_rme(const basis::RelativeStateLSJT& bra,
                                     const basis::RelativeStateLSJT& ket,
                                     const double osc_b) = 0;
        // Data Members
        Order order;
        int G0; // Parity
        int J0; // Tensor rank
        int T0; // Isotensor rank
    };
}

#endif
