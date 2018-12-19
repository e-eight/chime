#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_integration.h>
#include "constants.h"
#include "tdho.h"
#include "utility.h"
#include "operators.h"

/* 
   The struct data is used to pass the integration parameters to the GSL 
   integration routines.
*/
struct data { int n; int l; }; 

double radial_integrand(double y, void *userdata)
{
    /* 
       This is momentum radial integrand that appears in the calculation of the
       matrix elements of the D term. Essentially it appears twice, once with
       the ket parameters and once with the bra parameters.
       
       The integrand is exp(-y/2) y^((l+1)/2) L_n^(l+1/2)(y) and the domain of
       integration is [0, ∞). The integral is computed using Gauss-Laguerre
       quadrature, for which exp(-y/2) y^((l+1)/2) acts as the weight. Hence
       this function just needs to return the Laguerre polynomial part.
    */
    struct data *d = (struct data *) userdata;
    int n = d->n, l = d->l;
    return gsl_sf_laguerre_n(n, l + 0.5, y);
}

double radial_integral(int n, int l)
{
    /* 
       This integrates the above integral over the domain [0, ∞) using
       Gauss-Laguerre quadrature from the GSL.
    */

    const double alpha = (l + 1.0) / 2, beta = 0, gsl_a = 0, gsl_b = 0.5;
    const size_t num_nodes = n / 2 + 1;

    const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;
    gsl_integration_fixed_workspace *w
	= gsl_integration_fixed_alloc(T, num_nodes, gsl_a, gsl_b, alpha, beta);

    struct data params = { n, l };
    gsl_function g_func;
    g_func.function = &radial_integrand;
    g_func.params = &params;

    double result;
    gsl_integration_fixed(&g_func, &result, w);
    gsl_integration_fixed_free(w);

    return result;
}

double gamow_teller_n2lo_d_term(q_nums *ket, q_nums *bra, double *b)
{
    /* 
       This function computes the N2LO term proportional to the low energy
       constant D.
       
       This term contains two sinc functions, sinc(πl) and sinc(πl'), however
       the GSL implementation of sinc already assumes that π is present.
     */
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (abs(mtf - mti) != 1 || mjf != mji)
	return 0;
    
    double isospin_one = (2 * delta(tf, 1) * delta(ti, 1)
			  + delta(tf, 0) * delta(ti, 1)
			  - delta(tf, 1) * delta(ti, 0)); 

    double isospin_two = (2 * delta(tf, 1) * delta(ti, 1)
			  - delta(tf, 0) * delta(ti, 1)
			  + delta(tf, 1) * delta(ti, 0));

    double spin_one = (delta(sf, 1) * delta(si, 1)
		       + delta(sf, 0) * delta(si, 1)
		       + delta(sf, 1) * delta(si, 0));

    double spin_two = (delta(sf, 1) * delta(si, 1)
		       - delta(sf, 0) * delta(si, 1)
		       - delta(sf, 1) * delta(si, 0));

    double radial_term = radial_integral(ni, li) * radial_integral(nf, lf);
    radial_term *= (gsl_sf_sinc(lf) * gsl_sf_sinc(li)
		    * sqrt(2*li + 1) * sqrt(2*lf + 1));
    radial_term *= sqrt(exp(gsl_sf_lngamma(ni + 1)
			    + gsl_sf_lngamma(nf + 1)
			    - gsl_sf_lngamma(ni + li + 1.5)
			    - gsl_sf_lngamma(nf + lf + 1.5)));
    radial_term /= (4 * gsl_pow_2(M_PI) * gsl_pow_3(*b));

    double result = ((isospin_one * spin_one + isospin_two * spin_two)
		     * radial_term);
    result *= (gsl_pow_int(-1, ni + nf) * LEC_D
	       * cg_coeff(li, si, ji, 0, mji, mji)
	       * cg_coeff(lf, sf, jf, 0, mjf, mjf) / 4);
    result /= (cg_coeff(ti, 1, tf, 1, -1, 0) * cg_coeff(ji, 1, jf, mji, 0, mjf));
    
    return result;
}

double gamow_teller_n2lo(q_nums *ket, q_nums *bra, double *b)
{
    return gamow_teller_n2lo_d_term(ket, bra, b);
}

double gamow_teller_bare(q_nums *ket, q_nums *bra, double *b)
{
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (abs(mtf - mti) != 1)
	return 0;
    
    double result = 0;
    if (tf == 1 && ti == 1)
    {
	if (sf == 1 && si == 1)
	{
#pragma omp parallel for
	    for (int ms = -si; ms <= si; ms++)
	    {
		double cgproduct = (cg_coeff(li, si, ji, mji-ms, ms, mji)
				    * cg_coeff(lf, sf, jf, mjf-ms, ms, mjf));
		result += cgproduct * delta(li, lf) * delta(mji, mjf);
		result *= 2 * sqrt(2) * G_A * delta(ni, nf);
	    }
	}
    }
    else if (tf == 0 && ti == 1)
    {
	if ((sf == 0 && si == 1) || (sf == 1 && si == 0))
	{
	    result = (sqrt(2) * G_A * cg_coeff(li, si, ji, mji, 0, mji)
		      * cg_coeff(lf, sf, jf, mjf, 0, mjf) * delta(li, lf)
		      * delta(mji, mjf) * delta(ni, nf)); 
	}
    }
    else if (tf == 1 && ti == 0)
    {
	if ((sf == 0 && si == 1) || (sf == 1 && si == 0))
	{
	    result = - (sqrt(2) * G_A * cg_coeff(li, si, ji, mji, 0, mji)
		      * cg_coeff(lf, sf, jf, mjf, 0, mjf) * delta(li, lf)
		      * delta(mji, mjf) * delta(ni, nf));
	}
    }

    result /= (cg_coeff(ti, 1, tf, 1, -1, 0) * cg_coeff(ji, 1, jf, mji, 0, mjf));

    return result;
}

double operator_gt(char *order, q_nums *ket, q_nums *bra, double *b)
{
    if (strcmp(order, "bare") == 0)
	return gamow_teller_bare(ket, bra, b);
    else if (strcmp(order, "n2lo") == 0)
	return gamow_teller_n2lo(ket, bra, b);
    else if (strcmp(order, "all") == 0)
	return (gamow_teller_bare(ket, bra, b)
		+ gamow_teller_n2lo(ket, bra, b));
}
