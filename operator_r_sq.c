#include <math.h>
#include <string.h>
#include "constants.h"
#include "tdho.h"
#include "utility.h"
#include "operators.h"

double p_one(int l)
{
    return (2*l - 1) * (2*l + 1);
}

double p_two(int l, int m)
{
    return (l - m) * (l + m);
}

double p_three(int l, int m)
{
    return sqrt((l - m -1) * (l - m) * (l + m - 1) * (l + m));
}

double integral_y(int li, int mli, int lf, int mlf)
{
    if (mlf != mli)
    {
        return 0;
    }
    else if (lf == li + 2)
    {
        double coeff = sqrt((2*li + 5) / (2*li + 1));
        return coeff * p_three(li + 2, mli) / p_one(li + 1);
    }
    else if (lf == li - 2)
    {
        double coeff = sqrt((2*li - 3) / (2*li + 1));
        return coeff * p_three(li, mli) / p_one(li);
    }
    else if (lf == li)
    {
        double term1 = p_two(li, mli) / p_one(li);
        double term2 = p_two(li + 1, mli) / p_one(li + 1);
        return term1 + term2;
    }

    return 0;
}

double integral_l(int ni, int li, int nf, int lf, double b)
{
    double result = 0;
    
    if (lf == li + 2)
    {
        result = (sqrt((ni + li + 2.5) * (ni + li + 1.5)) * delta(nf, ni)
                  + 2 * sqrt(ni * (ni + li + 1.5)) * delta(nf, ni - 1)
                  + sqrt(ni * (ni - 1)) * delta(nf, ni - 2)) ;
    }
    else if (lf == li - 2)
    {
        result = (sqrt((ni + li + 0.5) * (ni + li - 0.5)) * delta(nf, ni)
                  + 2 * sqrt((ni + 1) * (ni + li - 0.5)) * delta(nf, ni - 1)
                  + sqrt((ni + 1) * (ni + 2)) * delta(nf, ni - 2));
    }
    else if (lf == li)
    {
        result = ((2*ni + li + 1.5) * delta(nf, ni)
                  + sqrt(ni * (ni + li + 0.5)) * delta(nf, ni - 1)
                  + sqrt((ni + 1) * (ni + li + 1.5)) * delta(ni, nf - 1));
    }

    return 6 * b * b * result;

}
    
double r_sq_lo(q_nums *ket, q_nums *bra, double *b)
{   
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (sf != si || tf != ti || mjf != mji || mtf != mti)
        return 0;

    double result = 0;
#pragma omp parallel for
    for (int ms = -si; ms <= si; ms++)
    {
        double cgproduct = (cg_coeff(li, si, ji, mji-ms, ms, mji)
                            * cg_coeff(lf, sf, jf, mjf-ms, ms, mjf));
        result += (cgproduct * integral_l(ni, li, nf, lf, *b)
                   * integral_y(li, mji-ms, lf, mjf-ms));
    }
    return 0.5 * result;
}

double r_sq_struct(q_nums *ket, q_nums *bra, double *b)
{
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (sf != si || tf != ti || mjf != mji || mtf != mti)
        return 0;
    
    double result = 0;
#pragma omp parallel for
    for (int ms = -si; ms <= si; ms++)
    {
        double cgproduct = (cg_coeff(li, si, ji, mji-ms, ms, mji)
                            * cg_coeff(lf, sf, jf, mjf-ms, ms, mjf));
        result += (cgproduct * delta(ni, nf) * delta(li, lf)
                   * delta(mji-ms, mjf-ms));
    }
    return 2 * R_ES_SQUARED * result;
}

double operator_r_sq(char *order, q_nums *ket, q_nums *bra, double *b)
{
    if (!strcmp(order, "LO"))
        return r_sq_lo(ket, bra, b);
    else if (!strcmp(order, "N2LO Struct"))
        return r_sq_struct(ket, bra, b);
}
