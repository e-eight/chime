#include <cmath>
#include <complex>
#include <gsl/gsl_sf_laguerre.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "wavefunction.h"

namespace tdho
{
    WaveFunction::WaveFunction():
	n(0), l(0), m(0), osc_b(0) {}
    WaveFunction::WaveFunction(double osc_b):
	n(0), l(0), m(0), osc_b(osc_b) {}
    WaveFunction::WaveFunction(int n, int l, int m, double osc_b):
	n(n), l(l), m(m), osc_b(osc_b) {}

    WaveFunction::~WaveFunction() {}

    double WaveFunction::norm_nl()
    {
	auto result = 0.5 * (std::log(2) - 3 * std::log(osc_b)
			     + boost::math::lgamma(n + 1)
			     - boost::math::lgamma(n + l + 1.5));
	return std::exp(result);
    }

    double WaveFunction::radial_nl(double r)
    {
	auto rb = r / osc_b;
	return (std::exp(- rb * rb / 2) * std::pow(rb, l)
		* gsl_sf_laguerre_n(n, l + 0.5, rb * rb));
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
