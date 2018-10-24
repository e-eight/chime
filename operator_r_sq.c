#include <math.h>
#include <string.h>
#include "constants.h"
#include "tdho.h"
#include "utility.h"
#include "operators.h"

double plm(int l, int m)
{
    if (l != 0 && l != m && l != -m)
        return sqrt((l - m) * (l + m) / (2*l + 1.0) / (2*l - 1.0));
    return 0;
}

double spherical_integral_lo(int li, int mli, int lf, int mlf)
{
    double result = 0;
    if (mlf == mli)
    {
        result = ((plm(li, mli) * plm(li, mli)
                   + plm(li + 1, mli) * plm(li + 1, mli)) * delta(lf, li)
                  + plm(li + 1, mli) * plm(li + 2, mli) * delta(lf, li + 2)
                  + plm(lf + 1, mli) * plm(lf + 2, mli) * delta(li, lf + 2));
    }

    return result;
}

double radial_integral_bare(int ni, int li, int nf, int lf, double b)
{
    double result = 0;
    if (lf == li)
    {
        result = ((2*ni + li + 1.5) * delta(nf, ni)
                  - sqrt((ni + 1) * (ni + li + 1.5)) * delta(nf, ni + 1)
                  - sqrt((nf + 1) * (nf + li + 1.5)) * delta(ni, nf + 1));
        result = b * b * result;
    }
    return result;
}

double radial_integral_lo(int ni, int li, int nf, int lf, double b)
{
    double result = 0;
    if (lf == li)
    {
        result = radial_integral_bare(ni, li, nf, lf, b);
    }
    else if (lf == li + 2)
    {
        result = (sqrt((nf + li + 2.5) * (nf + li + 3.5)) * delta(ni, nf)
                  + 2 * sqrt((nf + 1) / (nf + li + 2.5)) * delta(ni, nf + 1)
                  + (sqrt((nf + 1) * (nf + 2))
                     / (nf + li + 3.5) / (nf + li + 2.5)) * delta(ni, nf + 2));
        result = b * b * result;
    }
    else if (li == lf + 2)
    {
        result = radial_integral_lo(nf, lf, ni, li, b);
    }
    return result;
}

double r_sq_bare(q_nums *ket, q_nums *bra, double *b)
{
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (sf != si || jf != ji || tf != ti || mjf != mji || mtf != mti)
        return 0;

    double result = 0;
#pragma omp parallel for
    for (int ms = -si; ms <= si; ms++)
    {
        double cgproduct = (cg_coeff(li, si, ji, mji-ms, ms, mji)
                            * cg_coeff(lf, sf, jf, mjf-ms, ms, mjf));
        result += (cgproduct * delta(lf, li)
                   * radial_integral_bare(ni, li, nf, lf, *b));
    }
    return result;
}
    
double r_sq_lo(q_nums *ket, q_nums *bra, double *b)
{   
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (sf != si || jf != ji || tf != ti || mjf != mji || mtf != mti)
        return 0;

    double result = 0;
#pragma omp parallel for
    for (int ms = -si; ms <= si; ms++)
    {
        double cgproduct = (cg_coeff(li, si, ji, mji-ms, ms, mji)
                            * cg_coeff(lf, sf, jf, mjf-ms, ms, mjf));
        result += (cgproduct * radial_integral_lo(ni, li, nf, lf, *b)
                   * spherical_integral_lo(li, mji-ms, lf, mjf-ms));
    }
    return 0.75 * result;
}

double r_sq_n2lo(q_nums *ket, q_nums *bra, double *b)
{
    int ni = ket->n, li = ket->l, si = ket->s, ji = ket->j, ti = ket->t;
    int nf = bra->n, lf = bra->l, sf = bra->s, jf = bra->j, tf = bra->t;
    int mji = ket->mj, mjf = bra->mj, mti = ket->mt, mtf = bra->mt;

    if (sf != si || jf != ji || tf != ti || mjf != mji || mtf != mti)
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
    return R_ES_SQUARED * result;
}

double operator_r_sq(char *order, q_nums *ket, q_nums *bra, double *b)
{
    if (strcmp(order, "bare") == 0)
        return r_sq_bare(ket, bra, b);
    else if (strcmp(order, "lo") == 0)
        return r_sq_lo(ket, bra, b);
    else if (strcmp(order, "nlo") == 0)
        return 0;
    else if (strcmp(order, "n2lo") == 0)
        return r_sq_n2lo(ket, bra, b);
    else if (strcmp(order, "all") == 0)
        return r_sq_lo(ket, bra, b) + r_sq_n2lo(ket, bra, b);
}
