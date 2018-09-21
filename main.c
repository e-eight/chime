#include <stdio.h>
#include <stdlib.h>
#include "tdho.h"
#include "operator_id.h"

#define HEADING "ISU Two-body matrix elements\n"
#define OPERATOR "Operator: "
#define CONSTANTS "Relevant constants: "
#define QNUMS "Quantum numbers: "

int main() {
    int nmax = 200, lmax = 2, NMAX = 10, LMAX = 0;
    int s = 1, j = 1, mj = 0, t = 0, mt = 0;
    int sp = 1, jp = 1, mjp = 0, tp = 0, mtp = 0;

    double b = 1, result;
    q_nums *i_nums = malloc(sizeof(q_nums));
    q_nums *f_nums = malloc(sizeof(q_nums));

    printf(HEADING);
    printf(OPERATOR);
    printf("Toy\n");
    printf(CONSTANTS);
    printf("None\n");
    printf(QNUMS);
    printf("s' = %d, j' = %d, t' = %d, s = %d, j = %d, t = %d\n",
	   sp, jp, tp, s, j, t);
    printf("nmax = %d, lmax = %d, NMAX = %d, LMAX = %d", nmax, lmax, NMAX, LMAX);
    printf("\n");
    printf("N\tL\tn\tl\tNp\tLp\tnp\tlp\tMatrix Element\n");
    
    for (int L = 0, Lp = 0; L <= LMAX && Lp <= LMAX; L++, Lp++) {
	for (int l = 0, lp = 0; l <= lmax && lp <= lmax; l++, lp++) {
	    for (int N = 0, Np = 0; N <= NMAX && Np <= NMAX; N++, Np++) {
		for (int n = 0, np = 0; n <= nmax && np <= nmax; n++, np++) {
		    set_q_nums(i_nums, n, l, N, L, s, j, mj, t, mt);
		    set_q_nums(f_nums, np, lp, Np, Lp, sp, jp, mjp, tp, mtp);
		    result = operator_id(i_nums, f_nums, &b);
		    if (result)
			printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n",
			       N, L, n, l, Np, Lp, np, lp, result);
		}
	    }
	}
    }

    printf("-99999\n");
    free(f_nums);
    free(i_nums);
}
