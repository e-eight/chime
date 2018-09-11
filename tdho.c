#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include "tdho.h"

bool isequal_wf_params(wf_params i_params, wf_params f_params) {
    return i_params.n == f_params.n && i_params.l == f_params.l
	&& i_params.ml == f_params.ml && i_params.b == f_params.b;
}

bool isequal_q_nums(q_nums i_nums, q_nums f_nums) {
    int n = i_nums.n, N = i_nums.N, np = f_nums.n, Np = f_nums.N;
    int l = i_nums.l, L = i_nums.L, lp = f_nums.l, Lp = f_nums.L;
    int s = i_nums.s, j = i_nums.j, sp = f_nums.s, jp = f_nums.j;
    int t = i_nums.t, tp = f_nums.t;
    int mj = i_nums.mj, mt = i_nums.mt, mjp = f_nums.mj, mtp = f_nums.mt;

    return N == Np && n == np && l == lp && L == Lp && s == sp && j == jp
	&& mj == mjp && t ==tp && mt == mtp;
}

double radial_nl(double p, const wf_params* params) {
    int n = params->n;
    int l = params->l;
    double a = l + 0.5;

    double b = params->b;;
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

double norm_nl(const wf_params* params) {
    int n = params->n;
    int l = params->l;
    double b = params->b, result;

    result = log(2) + 3 * log(b)
	+ gsl_sf_lngamma(n+1) - gsl_sf_lngamma(n+l+1.5);
    
    result *= 0.5;
    
    return exp(result); 
}

double y_lm_real(double theta, double phi, const wf_params* params) {
    int l = params->l;
    int ml = params->ml;
    
    return gsl_sf_legendre_sphPlm(l, 0, cos(theta)) * cos(ml*phi);
}

double y_lm_imag(double theta, double phi, const wf_params* params) {
    int l = params->l;
    int ml = params->ml;
    
    return gsl_sf_legendre_sphPlm(l, 0, cos(theta)) * sin(ml*phi);
}

double wavefunc_nlm_real(double p, double theta,
			 double phi, const wf_params* params) {
    return norm_nl(params) * radial_nl(p, params)
	* y_lm_real(theta, phi, params);
}

double wavefunc_nlm_imag(double p, double theta,
			 double phi, const wf_params* params) {
    return norm_nl(params) * radial_nl(p, params)
	* y_lm_imag(theta, phi, params);
}
