#include "chiral.h"
#include "operator_rc_sq.h"

namespace chiral
{
    void calculate_rme(const ChiralOperator& op,
                       const basis::RelativeStateLSJT& bra,
                       const basis::RelativeStateLSJT& ket,
                       const double& osc_b,
                       double& rme)	
    {
        switch(op.name)
        {
        case ChiralOperator::Name::charge_radius:
            rme = rc_sq_rme(op, bra, ket, osc_b);
            break;
        // case ChiralOperator::Name::gamow_teller:
        //     rme = gt_rme(op, bra, ket, osc_b);
        //     break;
        default:
            rme = 1.0;
            break;
        }
    }
}
