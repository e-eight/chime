#include <stdio.h>
#include "wavefunction.h"

int main() {
    wf_params iparams = { 100, 0, 0, 1 };
    double N_nl;

    N_nl = norm_nl(iparams);

    printf("%f\n", N_nl);
    return 0;
}
