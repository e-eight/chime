#ifndef OPERATOR_RC_SQ_H
#define OPERATOR_RC_SQ_H

#include <string>
#include "chiral.h"
#include "lib/basis/lsjt_scheme.h"

namespace chiral
{
    class ChargeRadiusOperator: public Operator
    {
    private:
        double lo_matrix_element(const basis::RelativeStateLSJT,
                                 const basis::RelativeStateLSJT,
                                 const double);
        double nlo_matrix_element(const basis::RelativeStateLSJT,
                                  const basis::RelativeStateLSJT,
                                  const double);
        double n2lo_matrix_element(const basis::RelativeStateLSJT,
                                   const basis::RelativeStateLSJT,
                                   const double);
        double n3lo_matrix_element(const basis::RelativeStateLSJT,
                                   const basis::RelativeStateLSJT,
                                   const double);
    public:
        // Constructor
        ChargeRadiusOperator(): Operator("rc_sq") {}
        void set_matrix_element(const Order,
                                const basis::RelativeStateLSJT,
                                const basis::RelativeStateLSJT,
                                const double) override;
    };
}

#endif
