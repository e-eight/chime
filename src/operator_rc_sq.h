#ifndef OPERATOR_RC_SQ_H
#define OPERATOR_RC_SQ_H

#include "basis/lsjt_scheme.h"

namespace chiral
{
    struct ChiralOperator;
    
    double rc_sq_rme(const ChiralOperator& op,
                     const basis::RelativeStateLSJT& bra,
                     const basis::RelativeStateLSJT& ket,
                     const double& osc_b);
    
    double rc_sq_lo(const basis::RelativeStateLSJT& bra,
                    const basis::RelativeStateLSJT& ket,
                    const double& osc_b);
    
    double rc_sq_nlo(const basis::RelativeStateLSJT& bra,
                     const basis::RelativeStateLSJT& ket,
                     const double& osc_b);
    
    double rc_sq_n2lo(const basis::RelativeStateLSJT& bra,
                      const basis::RelativeStateLSJT& ket,
                      const double& osc_b);
    
    double rc_sq_n3lo(const basis::RelativeStateLSJT& bra,
                      const basis::RelativeStateLSJT& ket,
                      const double& osc_b);
    
    double rc_sq_n4lo(const basis::RelativeStateLSJT& bra,
                      const basis::RelativeStateLSJT& ket,
                      const double& osc_b);
    
    double rc_sq_full(const basis::RelativeStateLSJT& bra,
                      const basis::RelativeStateLSJT& ket,
                      const double& osc_b);
}

#endif
