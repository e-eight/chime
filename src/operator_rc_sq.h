#ifndef OPERATOR_RC_SQ_H
#define OPERATOR_RC_SQ_H

#include "basis/lsjt_scheme.h"
#include "chiral.h"

namespace chiral
{
    class ChargeRadiusOperator : public ChiralOperator
    {
    public:
        ChargeRadiusOperator();

        ~ChargeRadiusOperator();
        
        void calculate_rme(const basis::RelativeStateLSJT& bra,
                           const basis::RelativeStateLSJT& ket,
                           const double& osc_b,
                           double& rme) override;
    };

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
