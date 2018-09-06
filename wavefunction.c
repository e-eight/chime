#include <stdio.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <math.h>
#include "wavefunction.h"

double radial_nl(double p, const wf_params params) {
    int n = params.n;
    int l = params.l;
    double a = l + 0.5;

    double b = params.b;;
    double bp = b * p;
    double bp2 = bp * bp;

    double gaussian, power, laguerre_na;

    gaussian = exp(-bp2 / 2);

    power = gsl_pow_int(bp, l);	

    switch (n) {
    case '1':
	laguerre_na = gsl_sf_laguerre_1(a, bp2);
	break;
    case '2':
	laguerre_na = gsl_sf_laguerre_2(a, bp2);
	break;
    case '3':
	laguerre_na = gsl_sf_laguerre_3(a, bp2);
	break;
    default:
	laguerre_na = gsl_sf_laguerre_n(n, a, bp2);
	break;
    }

    return gaussian * power * laguerre_na; 
}

double norm_nl(const wf_params params) {
    int n = params.n;
    int l = params.l;
    double b = params.b, result;

    result = log(2) + 3 * log(b)
	+ gsl_sf_lngamma(n+1) - gsl_sf_lngamma(n+l+1.5);
    
    result *= 0.5;
    
    return exp(result); 
}

double y_lm_real(double theta, double phi, const wf_params params) {
    int l = params.l;
    int m = params.m;
    double ctheta = cos(theta);
    
    return gsl_sf_legendre_sphPlm(l, 0, ctheta) * cos(m*phi);
}

double y_lm_imag(double theta, double phi, const wf_params params) {
    int l = params.l;
    int m = params.m;
    double ctheta = cos(theta);
    
    return gsl_sf_legendre_sphPlm(l, 0, ctheta) * sin(m*phi);
}

double wavefunc_nlm_real(double p, double theta,
			 double phi, const wf_params params) {
    return norm_nl(params) * radial_nl(p, params)
	* y_lm_real(theta, phi, params);
}

double wavefunc_nlm_imag(double p, double theta,
			 double phi, const wf_params params) {
    return norm_nl(params) * radial_nl(p, params)
	* y_lm_imag(theta, phi, params);
}
