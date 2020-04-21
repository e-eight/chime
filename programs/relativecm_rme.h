#ifndef RELATIVECM_RME_H_
#define RELATIVECM_RME_H_

#include <array>
#include "basis/lsjt_operator.h"

namespace relcm {

  void ConstructMu2nNLOOperator(const basis::RelativeCMOperatorParametersLSJT& op_params,
                                const basis::RelativeCMSpaceLSJT& relcm_space,
                                std::array<basis::RelativeCMSectorsLSJT, 3>& relcm_sectors,
                                std::array<basis::OperatorBlocks<double>, 3>& relcm_matrices,
                                const double& oscillator_energy,
                                const double& R);

}

#endif
