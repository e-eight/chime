#include <cmath>
#include <algorithm>
// #include "basis/lsjt_scheme.h"
#include "constants.h"
#include "utility.h"
// #include "chiral.h"
#include "operator_rc_sq.h"

namespace chiral
{
    ChargeRadiusOperator::ChargeRadiusOperator():
	ChiralOperator(Name::charge_radius, 0, 0) {}

    ChargeRadiusOperator::~ChargeRadiusOperator() {}

    void ChargeRadiusOperator::calculate_rme(const basis::RelativeStateLSJT& bra,
					     const basis::RelativeStateLSJT& ket,
					     const double& osc_b,
					     double& rme)
    {
	switch(order)
	{
	case Order::lo:
	    rme = rc_sq_lo(bra, ket, osc_b);
	case Order::nlo:
	    rme = rc_sq_nlo(bra, ket, osc_b);
	case Order::n2lo:
	    rme = rc_sq_n2lo(bra, ket, osc_b);
	case Order::n3lo:
	    rme = rc_sq_n3lo(bra, ket, osc_b);
	case Order::n4lo:
	    rme = rc_sq_n4lo(bra, ket, osc_b);
	default:
	    rme = rc_sq_full(bra, ket, osc_b);
	}
    }
    
    double rc_sq_lo(const basis::RelativeStateLSJT& bra,
		    const basis::RelativeStateLSJT& ket,
		    const double& osc_b)
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
    
	double integral_nl = osc_b * osc_b;
	if (ni == nf)
	{
	    integral_nl *= 2 * ni + li + 1.5;
	}
	else if (std::abs(ni - nf) == 1)
	{
	    int n = std::min(ni, nf);
	    integral_nl *= - std::sqrt((n + 1) * (n + li + 1.5));
	}

	double clebsch_product = 0;
#pragma omp parallel for
	for (int ms = -si; ms <= si; ms++)
	{
	    clebsch_product += (util::clebsch(li, si, ji, mji-ms, ms, mji)
				* util::clebsch(lf, sf, jf, mjf-ms, ms, mjf));
	}
	
	return clebsch_product * integral_nl;		    
    }

    double rc_sq_nlo(const basis::RelativeStateLSJT& bra,
		     const basis::RelativeStateLSJT& ket,
		     const double& osc_b)
    {
	return 0;
    }

    double rc_sq_n2lo(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double& osc_b)
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

	double clebsch_product = 0;
#pragma omp parallel for
	for (int ms = -si; ms <= si; ms++)
	{
	    clebsch_product += (util::clebsch(li, si, ji, mji-ms, ms, mji)
				* util::clebsch(lf, sf, jf, mjf-ms, ms, mjf));
	}

	return constants::R_ES_SQUARED * clebsch_product;
    }

    double rc_sq_n3lo(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double& osc_b)
    {
	return 0;
    }

    double rc_sq_n4lo(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double& osc_b)
    {
	return 0;
    }

    double rc_sq_full(const basis::RelativeStateLSJT& bra,
		      const basis::RelativeStateLSJT& ket,
		      const double& osc_b)
    {
	return (rc_sq_lo(bra, ket, osc_b)
		+ rc_sq_nlo(bra, ket, osc_b)
		+ rc_sq_n2lo(bra, ket, osc_b)
		+ rc_sq_n3lo(bra, ket, osc_b)
		+ rc_sq_n4lo(bra, ket, osc_b));
    }
    
}
