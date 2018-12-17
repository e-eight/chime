#ifndef CHARGE_RADIUS_H
#define CHARGE_RADIUS_H

#include <string>
#include "chiral.h"
#include "lib/basis/lsjt_scheme.h"

namespace chiral
{
    class ChargeRadiusOperator: public Operator
    {
    public:
        // Constructor
        ChargeRadiusOperator(): Operator("rc_sq") {}
        void set_matrix_element(const Order,
                                const basis::RelativeStateLSJT,
                                const basis::RelativeStateLSJT,
                                const double);
    };
}

#endif
