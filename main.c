#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "utility.h"
#include "tdho.h"
#include "operators.h"

int main(int argc, char *argv[])
{
    int nmax = 7, lmax = 2, NMAX = 0, LMAX = 0;
    int si = 1, ji = 1, mji = 0, ti = 0, mti = 0;
    int sf = 1, jf = 1, mjf = 0, tf = 0, mtf = 0;

    double osc_energy, osc_constant;

    osc_energy = 20; // in MeV
    osc_constant = HBARC / sqrt(RED_NUCLEON_MASS * osc_energy);
    
    q_nums *ket = malloc(sizeof(q_nums));
    q_nums *bra = malloc(sizeof(q_nums));

    printf("nmax = %d, lmax = %d \n", 2*nmax + lmax, lmax);
    printf(" ki  li  kf  lf   RME\n");
    
    double result;
    // for (int L = 0, Lp = 0; L <= LMAX && Lp <= LMAX; L++, Lp++) {
    for (int li = 0; li <= lmax; li = li + 2)
    {
        for (int lf = 0; lf <= lmax; lf = lf + 2)
        {
            //	    for (int N = 0, Np = 0; N <= NMAX && Np <= NMAX; N++, Np++) {
            for (int ni = 0; ni <= nmax; ni++)
            {
                for (int nf = 0; nf <= nmax; nf++)
                {
                    set_q_nums(ket, ni, li, 0, 0, si, ji, mji, ti, mti);
                    set_q_nums(bra, nf, lf, 0, 0, sf, jf, mjf, tf, mtf);
                    result = operator_r_sq("N2LO Struct", ket, bra, &osc_constant);
                    printf("%3d %3d %3d %3d   %e\n",
                           2*ni + li, li, 2*nf + lf, lf, result);
                }
            }
        }
    }

    printf("-99999\n");
    free(bra);
    free(ket);
}
