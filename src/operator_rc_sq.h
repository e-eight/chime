#ifndef OPERATOR_RC_SQ_H
#define OPERATOR_RC_SQ_H

#include "chiral.h"
#include "basis/lsjt_scheme.h"

namespace chiral
{
    double operator_rc_sq(Operator::Orders order,
                          basis::RelativeStateLSJT bra,
                          basis::RelativeStateLSJT ket,
                          double osc_b);
    
    double rc_sq_lo(basis::RelativeStateLSJT bra,
                    basis::RelativeStateLSJT ket,
                    double osc_b);
    
    double rc_sq_nlo(basis::RelativeStateLSJT bra,
                     basis::RelativeStateLSJT ket,
                     double osc_b);
    
    double rc_sq_n2lo(basis::RelativeStateLSJT bra,
                      basis::RelativeStateLSJT ket,
                      double osc_b);
    
    double rc_sq_n3lo(basis::RelativeStateLSJT bra,
                      basis::RelativeStateLSJT ket,
                      double osc_b);
    
    double rc_sq_n4lo(basis::RelativeStateLSJT bra,
                      basis::RelativeStateLSJT ket,
                      double osc_b);
    
    double rc_sq_full(basis::RelativeStateLSJT bra,
                      basis::RelativeStateLSJT ket,
                      double osc_b);
}

#endif
