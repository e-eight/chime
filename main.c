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
    

    char *operator, *order, *osc_energy, *nmax;
    char *lmin_i, *lmax_i, *lmin_f, *lmax_f, *ji, *jf;
    double osc_constant;

    operator = argv[1];
    order = argv[2];
    osc_energy = argv[3];
    nmax = argv[4];
    lmin_i = argv[5];
    lmax_i = argv[6];
    lmin_f = argv[7];
    lmax_f = argv[8];
    ji = argv[9];
    jf = argv[10];

    int si = 1, mji = 0, ti = 0, mti = 0;
    int sf = 1, mjf = 0, tf = 0, mtf = 0;

    osc_constant = HBARC / sqrt(RED_NUCLEON_MASS * atof(osc_energy));
    
    q_nums *ket = malloc(sizeof(q_nums));
    q_nums *bra = malloc(sizeof(q_nums));

    // printf("nmax = %d \n", 2 * atoi(nmax) + atoi(lmax));
    printf(" ni  li  nf  lf   RME\n");
    
    double result;
    // for (int L = 0, Lp = 0; L <= LMAX && Lp <= LMAX; L++, Lp++) {
    for (int li = atoi(lmin_i); li <= atoi(lmax_i); li = li + 2)
    {
        for (int lf = atoi(lmin_f); lf <= atoi(lmax_f); lf = lf + 2)
        {
            //	    for (int N = 0, Np = 0; N <= NMAX && Np <= NMAX; N++, Np++) {
            for (int ni = 0; ni <= atoi(nmax); ni++)
            {
                for (int nf = 0; nf <= atoi(nmax); nf++)
                {
                    set_q_nums(ket, ni, li, 0, 0, si, atoi(ji), mji, ti, mti);
                    set_q_nums(bra, nf, lf, 0, 0, sf, atoi(jf), mjf, tf, mtf);

                    if (!strcmp(operator, "rsq"))
                    {
                        result = operator_r_sq(order, ket, bra, &osc_constant);
                        printf("%3d %3d %3d %3d   %e\n",
                               2*ni + li, li, 2*nf + lf, lf, result);
                    }
                }
            }
        }
    }

    printf("-99999\n");
    free(bra);
    free(ket);
}
