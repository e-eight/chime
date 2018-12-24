#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include "basis/lsjt_scheme.h"
#include "constants.h"
#include "utility.h"
#include "chiral.h"
#include "operator_rc_sq.h"
#include "wavefunction.h"

namespace bmq = boost::math::quadrature;
namespace chiral
{
    ChargeRadiusOperator::ChargeRadiusOperator():
	ChiralOperator() {}

    ChargeRadiusOperator::~ChargeRadiusOperator() {}

    double
    ChargeRadiusOperator::calculate_rme(const basis::RelativeStateLSJT& bra,
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
    
/********************************************************************************
                                 LO contribution
********************************************************************************/
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
                                NLO contribution
********************************************************************************/
    
    double rc_sq_nlo(const basis::RelativeStateLSJT& bra,
		     const basis::RelativeStateLSJT& ket,
		     const double osc_b)
    {
	return 0;
    }

/********************************************************************************
                               N2LO contribution
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
                               N3LO contribution
********************************************************************************/

    // Isospin matrix element of \vec{τ}_1 ⋅ \vec{τ}_2.
    double tau_dot_product(const int ti, const int tf)
    {
	if (ti != tf)
	    return 0;
	else
	    return (-3 * util::delta(ti, 0) + util::delta(ti, 1));
    }

    // Spin matrix element of \vec{σ}_1 ⋅ \hat{z} \vec{σ}_2 ⋅ \hat{z}.
    double sigma_z_dot_product(const int si, const int sf,
			       const int msi, const int msf)
    {
	if (si != sf && msi != msf)
	    return 0;
	else
	    return (-util::delta(msi, 0) + util::delta(msi, 1)
		    + util::delta(msi, -1));
    }

    double radial_integrand_n3lo(tdho::WaveFunction& wf1,
				 tdho::WaveFunction& wf2,
				 double r)
    {
	if (wf1.l != wf2.l)
	    return 0;
	
	auto norm1 = wf1.norm_nl(), norm2 = wf2.norm_nl();
	auto radial1 = wf1.radial_nl(r), radial2 = wf2.radial_nl(r);
	return (norm1 * norm2 * radial1 * radial2
		* std::exp(-constants::PION_MASS * r) * r);
    }
    
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

	bool kronecker = (lf == lf && sf == si && jf == ji && tf == ti
			  && mjf == mji && mtf == mti);
	if (!kronecker)
	    return 0;

	double result = 0;
	auto wf1 = tdho::WaveFunction(osc_b), wf2 = tdho::WaveFunction(osc_b);
	wf1.n = ni; wf1.l = li;
	wf2.n = nf; wf2.l = lf;
#pragma omp parallel for
	for (int ms = -si; ms <= si; ms++)
	{
	    auto cg_product = (util::clebsch(li, si, ji, mji - ms, ms, mji)
			       * util::clebsch(lf, sf, jf, mjf - ms, ms, mjf));
	    wf1.m = mji - ms; wf2.m = mjf - ms;
	    auto integrand =
		[&](double r) { return radial_integrand_n3lo(wf1, wf2, r); };
	    double error, Q;
	    Q = bmq::gauss_kronrod
		<double, 15>::integrate(integrand, 0,
					std::numeric_limits<double>::infinity(),
					5, 1e-9, &error);
	    result += cg_product * Q;
	}

	result *= (-0.75 * tau_dot_product(ti, tf)
		   * std::pow(constants::G_A / constants::F_PION, 2)
		   / (2 * constants::RED_NUCLEON_MASS));
	
	return result;
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
