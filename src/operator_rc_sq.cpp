#include <cmath>
#include <algorithm>
#include "constants.h"
#include "utility.h"
#include "basis/lsjt_scheme.h"
#include "operator_rc_sq.h"

double chiral::operator_rc_sq(chiral::Operator::Order order,
			      basis::RelativeStateLSJT bra,
			      basis::RelativeStateLSJT ket,
			      double osc_b)
{
    switch(order)
    {
    case chiral::Operator::Order::lo:
	rme = rc_sq_lo(bra, ket, osc_b);
	break;
    case chiral::Operator::Order::nlo:
	rme = rc_sq_nlo(bra, ket, osc_b);
	break;
    case chiral::Operator::Order::n2lo:
	rme = rc_sq_n2lo(bra, ket, osc_b);
	break;
    case chiral::Operator::Order::n3lo:
	rme = rc_sq_n3lo(bra, ket, osc_b);
	break;
    case chiral::Operator::Order::n4lo:
	rme = rc_sq_n4lo(bra, ket, osc_b);
	break;
    default:
	rme = rc_sq_full(bra, ket, osc_b);
	break;
    }
}

double chiral::rc_sq_lo(basis::RelativeStateLSJT bra,
			basis::RelativeStateLSJT ket,
			double osc_b)
{
    int ni = ket.N(), li = ket.L(), si = ket.S(), ji = ket.J(), ti = ket.T();
    int nf = bra.N(), lf = bra.L(), sf = bra.S(), jf = bra.J(), tf = bra.T();
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

    return clebsch_prouct * integral_nl;		    
}

double chiral::rc_sq_nlo(basis::RelativeStateLSJT bra,
			 basis::RelativeStateLSJT ket,
			 double osc_b)
{
    return 0;
}

double chiral::rc_sq_n2lo(basis::RelativeStateLSJT bra,
			  basis::RelativeStateLSJT ket,
			  double osc_b)
{
    int ni = ket.N(), li = ket.L(), si = ket.S(), ji = ket.J(), ti = ket.T();
    int nf = bra.N(), lf = bra.L(), sf = bra.S(), jf = bra.J(), tf = bra.T();
    int mji = 0, mjf = 0, mti = 0, mtf = 0;

    bool kronecker = (ni == nf && li == lf && sf == si && jf == ji && tf == ti
		      && mjf == mji && mtf == mti);
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

double chiral::rc_sq_n3lo(basis::RelativeStateLSJT bra,
			  basis::RelativeStateLSJT ket,
			  double osc_b)
{
    return 0;
}

double chiral::rc_sq_n4lo(basis::RelativeStateLSJT bra,
			  basis::RelativeStateLSJT ket,
			  double osc_b)
{
    return 0;
}

double chiral::rc_sq_full(basis::RelativeStateLSJT bra,
			  basis::RelativeStateLSJT ket,
			  double osc_b)
{
    return (chiral::rc_sq_lo(bra, ket, osc_b)
	    + chiral::rc_sq_nlo(bra, ket, osc_b)
	    + chiral::rc_sq_n2lo(bra, ket, osc_b)
	    + chiral::rc_sq_n3lo(bra, ket, osc_b)
	    + chiral::rc_sq_n4lo(bra, ket, osc_b));
}
