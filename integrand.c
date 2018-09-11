#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "tdho.h"
#include "cuba.h"

double radial_1d(double p, void* params) {
    return radial_nl(p, (wf_params*) params);
}
