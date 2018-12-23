#include <cmath>
#include "basis/lsjt_scheme.h"
#include "constants.h"
#include "utility.h"
#include "chiral.h"
#include "operator_rc_sq.h"
#include "wavefunction.h"

namespace chiral
{
    ChargeRadiusOperator::ChargeRadiusOperator():
	ChiralOperator() {}

    ChargeRadiusOperator::~ChargeRadiusOperator() {}

    double ChargeRadiusOperator::calculate_rme(const basis::RelativeStateLSJT& bra,
					       const basis::RelativeStateLSJT& ket,
					       const double osc_b)
    {
	switch(order)
	{
	case Order::lo:
	    return rc_sq_lo(bra, ket, osc_b);
	case Order::nlo:
	    return rc_sq_nlo(bra, ket, osc_b);
	case Order::n2lo:
	    return rc_sq_n2lo(bra, ket, osc_b);
	case Order::n3lo:
	    return rc_sq_n3lo(bra, ket, osc_b);
	case Order::n4lo:
	    return rc_sq_n4lo(bra, ket, osc_b);
	default:
	    return rc_sq_full(bra, ket, osc_b);
	}
    }
    
/*******************************************************************************
                         Leading order contribution
*******************************************************************************/
    double radial_integral_lo(int ni, int li, int nf, int lf, const double osc_b)
    {
	if (lf == li)
	{
	    auto result = ((2 * ni + li + 1.5) * util::delta(nf, ni)
			   - (std::sqrt((ni + 1) *  (ni + li + 1.5))
			      * util::delta(nf, ni + 1))
			   - (std::sqrt(nf + 1) * (nf + li + 1.5)
			      * util::delta(ni, nf + 1)));
	    return osc_b * osc_b * result;
	}
	else
	    return 0;
    }
 
    double rc_sq_lo(const basis::RelativeStateLSJT& bra,
		    const basis::RelativeStateLSJT& ket,
		    const double osc_b)
    {
	int ni = ket.N(), nf = bra.N();
	int li = ket.L(), lf = bra.L();
	int si = ket.S(), sf = bra.S();
	int ji = ket.J(), jf = bra.J();
	int ti = ket.T(), tf = bra.T();
	int mji = 0, mjf = 0, mti = 0, mtf = 0;

	bool kronecker = (li == lf && sf == si && jf == ji && tf == ti
		      && mjf == mji && mtf == mti);
	if (!kronecker)
	    return 0;
	else if (std::abs(ni - nf) > 1)
	    return 0;
    
	double result = 0;
#pragma omp parallel for
	for (int ms = -si; ms <= si; ms++)
	{
	    auto cg_product = (util::clebsch(li, si, ji, mji-ms, ms, mji)
			       * util::clebsch(lf, sf, jf, mjf-ms, ms, mjf));
	    result += cg_product * radial_integral_lo(ni, li, nf, lf, osc_b);
	}
	
	return result;		    
    }

/********************************************************************************
                        Next-to-leading order contribution
********************************************************************************/
    
    double rc_sq_nlo(const basis::RelativeStateLSJT& bra,
		     const basis::RelativeStateLSJT& ket,
		     const double osc_b)
    {
	return 0;
    }

/********************************************************************************
                   Next-to-next-to-leading order contribution
********************************************************************************/

    double rc_sq_n2lo(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double osc_b)
    {
	int ni = ket.N(), nf = bra.N();
	int li = ket.L(), lf = bra.L();
	int si = ket.S(), sf = bra.S();
	int ji = ket.J(), jf = bra.J();
	int ti = ket.T(), tf = bra.T();
	int mji = 0, mjf = 0, mti = 0, mtf = 0;
    
	bool kronecker = (ni == nf && li == lf && sf == si && jf == ji
			  && tf == ti && mjf == mji && mtf == mti);
	if (!kronecker)
	    return 0;

	double cg_product = 0;
#pragma omp parallel for
	for (int ms = -si; ms <= si; ms++)
	{
	    cg_product += (util::clebsch(li, si, ji, mji-ms, ms, mji)
			   * util::clebsch(lf, sf, jf, mjf-ms, ms, mjf));
	}

	return constants::R_ES_SQUARED * cg_product;
    }

/********************************************************************************
              Next-to-next-to-next-to-leading order contribution
********************************************************************************/

    double rc_sq_n3lo(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double osc_b)
    {
	int ni = ket.N(), nf = bra.N();
	int li = ket.L(), lf = bra.L();
	int si = ket.S(), sf = bra.S();
	int ji = ket.J(), jf = bra.J();
	int ti = ket.T(), tf = bra.T();
	int mji = 0, mjf = 0, mti = 0, mtf = 0;

	bool kronecker = (sf == si && jf == ji && tf == ti
			  && mjf == mji && mtf == mti);
	if (!kronecker)
	    return 0;
	
	return 0;
    }

    double rc_sq_n4lo(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double osc_b)
    {
	return 0;
    }

    double rc_sq_full(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double osc_b)
    {
	return (rc_sq_lo(bra, ket, osc_b)
		+ rc_sq_nlo(bra, ket, osc_b)
		+ rc_sq_n2lo(bra, ket, osc_b)
		+ rc_sq_n3lo(bra, ket, osc_b)
		+ rc_sq_n4lo(bra, ket, osc_b));
    }
    
}
