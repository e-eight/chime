#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include "tdho.h"

void set_q_nums(q_nums *nums, int n, int l, int N, int L,
		int s, int j, int mj, int t, int mt)
{
    nums->n = n;
    nums->l = l;
    nums->N = N;
    nums->L = L;
    nums->s = s;
    nums->j = j;
    nums->mj = mj;
    nums->t = t;
    nums->mt = mt;
}

bool isequal_wf_params(wf_params *ket, wf_params *bra)
{
    return (ket->n == bra->n && ket->l == bra->l
	    && ket->ml == bra->ml && ket->b == bra->b);
}

bool isequal_q_nums(q_nums *ket, q_nums *bra)
{
    int n = ket->n, N = ket->N, np = bra->n, Np = bra->N;
    int l = ket->l, L = ket->L, lp = bra->l, Lp = bra->L;
    int s = ket->s, j = ket->j, sp = bra->s, jp = bra->j;
    int t = ket->t, tp = bra->t;
    int mj = ket->mj, mt = ket->mt, mjp = bra->mj, mtp = bra->mt;

    return (N == Np && n == np && l == lp && L == Lp && s == sp && j == jp
	    && mj == mjp && t ==tp && mt == mtp);
}

double radial_nl(double p, const wf_params *params)
{
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

double norm_nl(const wf_params *params)
{
    int n = params->n;
    int l = params->l;
    double b = params->b, result;

    result = (log(2) + 3 * log(b)
	      + gsl_sf_lngamma(n+1) - gsl_sf_lngamma(n+l+1.5));
    
    result *= 0.5;
    
    return exp(result); 
}

double y_lm_real(double theta, double phi, const wf_params *params)
{
    int l = params->l;
    int ml = params->ml;
    
    return gsl_sf_legendre_sphPlm(l, 0, cos(theta)) * cos(ml*phi);
}

double y_lm_imag(double theta, double phi, const wf_params *params)
{
    int l = params->l;
    int ml = params->ml;
    
    return gsl_sf_legendre_sphPlm(l, 0, cos(theta)) * sin(ml*phi);
}

double wavefunc_nlm_real(double p, double theta,
			 double phi, const wf_params *params)
{
    return (norm_nl(params) * radial_nl(p, params)
	    * y_lm_real(theta, phi, params));
}

double wavefunc_nlm_imag(double p, double theta,
			 double phi, const wf_params *params)
{
    return (norm_nl(params) * radial_nl(p, params)
	    * y_lm_imag(theta, phi, params));
}
