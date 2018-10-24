#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "utility.h"
#include "tdho.h"
#include "operators.h"

void print_rme(char *operator, char *order, q_nums *ket,
               q_nums *bra, double *osc_constant)
{
    double result = 0;
    int ni = ket->n, li = ket->l, nf = bra->n, lf = bra->l;
    if (strcmp(operator, "rsq") == 0)
        result = operator_r_sq(order, ket, bra, osc_constant);
    printf("%3d %3d %3d %3d   %e\n", 2*ni + li, li, 2*nf + lf, lf, result);
}

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

    int si = 1, ti = 0, mji = 0, mti = 0;
    int sf = 1, tf = 0, mjf = 0, mtf = 0;

    osc_constant = HBARC / sqrt(RED_NUCLEON_MASS * atof(osc_energy));
    
    q_nums *ket = malloc(sizeof(q_nums));
    q_nums *bra = malloc(sizeof(q_nums));

    printf(" ni  li  nf  lf   RME\n");
    
    for (int li = atoi(lmin_i); li <= atoi(lmax_i); li = li + 2)
    {
        for (int lf = atoi(lmin_f); lf <= atoi(lmax_f); lf = lf + 2)
        {
            for (int ni = 0; ni <= atoi(nmax); ni = ni + 2)
            {
                for (int nf = 0; nf <= atoi(nmax); nf = nf + 2)
                {
                    if (ni >= li && nf >= lf)
                    {
                        int nri = (ni - li) / 2;
                        int nrf = (nf - lf) / 2;
                        set_q_nums(ket, nri, li, 0, 0, si,
                                   atoi(ji), mji, ti, mti);
                        set_q_nums(bra, nrf, lf, 0, 0, sf,
                                   atoi(jf), mjf, tf, mtf);
                        print_rme(operator, order, ket, bra, &osc_constant);
                    }
                }
            }
        }
    }

    printf("-99999\n");
    free(bra);
    free(ket);
}
