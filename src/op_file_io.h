#include "chiral.h"

namespace chiral
{
    struct MatrixParameters
    {
	ChiralOperator op;
	int J0_min = 0, J0_max = 0; // Limits of tensorial order of operator.
	int T0_min = 0, T0_max = 0; // Limits of isotensorial order of operator.
	int Nmax = 0, Jmax = 0; // Cutoff parameters for element generation.
	double hw = 0; // Oscillator energy.
    };

    void read_input_to_params(MatrixParameters&, const int, const int);
}
