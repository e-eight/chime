#ifndef RELATIVE_RME_H_
#define RELATIVE_RME_H_

#include <array>

#include "basis/lsjt_operator.h"

namespace chime {
namespace relative {

void ConstructMu2nNLOOperator(
    const basis::RelativeOperatorParametersLSJT& op_params,
    const basis::RelativeSpaceLSJT& rel_space,
    std::array<basis::RelativeSectorsLSJT, 3>& rel_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& rel_matrices,
    const double& oscillator_energy, const double& R);

}  // namespace relative
}  // namespace chime

#endif
