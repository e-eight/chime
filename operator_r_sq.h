#ifndef __OPERATOR_R_SQ_H__
#define __OPERATOR_R_SQ_H__

typedef struct integrand_l_params integrand_l_params;
struct integrand_l_params { int n; double alpha; int np; double alphap; };
double integral_l1(wf_params *i_params, wf_params *f_params);

#endif
