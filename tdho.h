#ifndef __TDHO_H__
#define __TDHO_H__

typedef struct wf_params wf_params;
struct wf_params { int n; int l; int ml; double b; };

typedef struct q_nums q_nums;
struct q_nums { int n; int l; int N; int L; int s;
    int j; int mj; int t; int mt; };

bool isequal_wf_params(wf_params i_params, wf_params f_params);
bool isequal_q_nums(q_nums i_nums, q_nums f_nums);
double radial_nl(double p, const wf_params* params);
double norm_nl(const wf_params* params);
double y_lm_real(double theta, double phi, const wf_params* params);
double y_lm_imag(double theta, double phi, const wf_params* params);
double wavefunc_nlm_real(double p, double theta,
			 double phi, const wf_params* params);
double wavefunc_nlm_imag(double p, double theta,
			 double phi, const wf_params* params);

#endif
