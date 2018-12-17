#ifndef OPERATOR_RC_SQ_H
#define OPERATOR_RC_SQ_H

#include "chiral.h"
#include "lib/basis/lsjt_scheme.h"

namespace chiral
{
    class ChargeRadiusOperator: public Operator
    {
    public:
        ChargeRadiusOperator(): Operator(Operator::Names::charge_radius) {}
        void set_rme() override;
    };

    double rc_sq_lo_rme(basis::RelativeStateLSJT bra,
                        basis::RelativeStateLSJT ket,
                        double osc_b);
    double rc_sq_nlo_rme(basis::RelativeStateLSJT bra,
                         basis::RelativeStateLSJT ket,
                         double osc_b);
    double rc_sq_n2lo_rme(basis::RelativeStateLSJT bra,
                          basis::RelativeStateLSJT ket,
                          double osc_b);
    double rc_sq_n3lo_rme(basis::RelativeStateLSJT bra,
                          basis::RelativeStateLSJT ket,
                          double osc_b);
    double rc_sq_n4lo_rme(basis::RelativeStateLSJT bra,
                          basis::RelativeStateLSJT ket,
                          double osc_b);
    double rc_sq_full_rme(basis::RelativeStateLSJT bra,
                          basis::RelativeStateLSJT ket,
                          double osc_b);
}

#endif
