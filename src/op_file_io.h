#include <tuple>
#include <string>
#include "chiral.h"
#include "operator_rc_sq.h"

namespace chiral
{
    struct MatrixParameters
    {
        std::string name;
        std::string order;
        ChiralOperator* op;
        int J0_min, J0_max; // Tensor limits
        int T0_min, T0_max; // Isotensor limits
        int Nmax, Jmax;     // Basis cutoff
        double hw;          // Oscillator energy
    };
    
    MatrixParameters input_to_params();
}
