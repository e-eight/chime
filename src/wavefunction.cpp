#include <cmath>
#include <complex>
#include <gsl/gsl_sf_laguerre.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "wavefunction.h"

namespace tdho
{
    WaveFunction::WaveFunction():
	osc_b(0), n(0), l(0), m(0) {}
    WaveFunction::WaveFunction(double osc_b):
	osc_b(osc_b), n(0), l(0), m(0) {}
    WaveFunction::WaveFunction(double osc_b, int n, int l, int m):
	osc_b(osc_b), n(n), l(l), m(m) {}

    WaveFunction::~WaveFunction() {}

    double WaveFunction::norm_nl()
    {
	auto result = 0.5 * (std::log(2) + 3 * std::log(3)
			     + boost::math::lgamma(n + 1)
			     - boost::math::lgamma(n + l + 1.5));
	return std::pow(-1, n) * std::exp(result);
    }

    double WaveFunction::radial_nl(double p)
    {
	auto bp = osc_b * p;
	return (std::exp(- bp * bp / 2) * std::pow(bp, l)
		* gsl_sf_laguerre_n(n, l + 0.5, bp * bp));
    }

    double WaveFunction::theta_lm(double theta)
    {
	// Theta part of the spherical harmonics.
	return boost::math::spherical_harmonic_r(l, m, theta, 0);
    }

    std::complex<double> WaveFunction::spherical_lm(double theta, double phi)
    {
	return boost::math::spherical_harmonic(l, m, theta, phi);
    }
}
