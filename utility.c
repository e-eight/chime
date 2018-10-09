#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_coupling.h>
#include "utility.h"

extern int delta(int a, int b);

double cg_coeff (int ja, int jb, int j, int mja, int mjb, int mj)
{
    int triangle_cond, mj_cond;
    triangle_cond = abs(ja-jb) <= j && j <= (ja+jb);
    mj_cond = mja + mjb == mj;
    if (!triangle_cond || !mj_cond)
	return 0;
    else
	return (gsl_pow_int(-1, ja-jb+mj) * sqrt(2*j + 1)
		* gsl_sf_coupling_3j(2*ja, 2*jb, 2*j, 2*mja, 2*mjb, -2*mj));
}

void print_preamble(char *my_operator, char *my_constants)
{
    printf(HEADING);
    printf(OPERATOR);
    printf(strcat(my_operator, "\n"));
    printf(CONSTANTS);
    printf(strcat(my_constants, "\n"));
    printf(QNUMS);
}
