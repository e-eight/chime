#include "relative_rme.h"

#include <Eigen/Dense>
#include <cmath>
#include <vector>

#include "basis_func/ho.h"
#include "chime.h"
#include "constants.h"
#include "quadpp/quadpp.h"
#include "quadpp/spline.h"
#include "tprme.h"

namespace chime {
namespace relative {

constexpr double mPi = constants::pion_mass_fm;
constexpr double mN = constants::nucleon_mass_fm;
constexpr double FPi = constants::pion_decay_constant_fm;
constexpr double gA = constants::gA;

///////////////////////////////////////////////////////////////////////////
/////////////////////// Magnetic moment matrix element ////////////////////
///////////////////////////////////////////////////////////////////////////

void ConstructMu2nNLOOperator(
    const basis::RelativeOperatorParametersLSJT& op_params,
    const basis::RelativeSpaceLSJT& rel_space,
    std::array<basis::RelativeSectorsLSJT, 3>& rel_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& rel_matrices,
    const double& oscillator_energy, const double& R)
{
  std::cout << " Constructing M1 operator...\n";
  assert(op_params.J0 == 1);
  assert(op_params.g0 == 0);
  assert((op_params.T0_min == 1) && (op_params.T0_max == 1));

  // Alias isospin rank.
  int T0 = op_params.T0_min;

  // Generate required harmonic oscillator basis functions, and radial
  // integral weights.
  int Nmax = op_params.Nmax;

  int npts = 3001;
  const double low = 0, high = 1;
  Eigen::ArrayXd x, r, jac, wt;
  quadpp::SemiInfiniteIntegralMesh(npts, low, high, x, r, jac);
  wt = r * r * jac;  // weights for radial integral with transformed variable

  std::cout << "  Generating basis functions...\n";
  std::vector<Eigen::ArrayXXd> ho_wfs;
  double brel = chime::RelativeOscillatorLength(oscillator_energy);
  basis_func::ho::WaveFunctionsUptoMaxL(ho_wfs, r, Nmax, Nmax, brel,
                                        basis_func::Space::coordinate);

  // Radial integral kernels.
  std::cout << "  Generating integral kernels...\n";
  Eigen::ArrayXd mpir = mPi * r;
  Eigen::ArrayXd expmpir = Eigen::exp(-mpir);
  Eigen::ArrayXd ypir = expmpir / mpir;
  Eigen::ArrayXd zpir = (1. + mpir);
  Eigen::ArrayXd tpir = (-1. + 2 * mpir);

  // Semilocal coordinate space regulator.
  Eigen::ArrayXd scs_reg =
      r.unaryExpr([&R](double rs) { return chime::SCSRegulator(rs, R); });

  // Zero initialize operator.
  std::cout << "  Zero initializing operator...\n";
  basis::ConstructZeroOperatorRelativeLSJT(op_params, rel_space, rel_sectors,
                                           rel_matrices);

  // Select T0 component.
  const basis::RelativeSectorsLSJT& sectors = rel_sectors[T0];
  basis::OperatorBlocks<double>& matrices = rel_matrices[T0];

  // Reduced matrix element calculation.
  std::cout << "  Starting matrix element calculation...\n";
  for (std::size_t sector_index = 0; sector_index < sectors.size();
       ++sector_index) {
    const basis::RelativeSectorsLSJT::SectorType& sector =
        sectors.GetSector(sector_index);
    const basis::RelativeSubspaceLSJT& bra_subspace = sector.bra_subspace();
    const basis::RelativeSubspaceLSJT& ket_subspace = sector.ket_subspace();

    // Alias for matrix.
    basis::OperatorBlock<double>& matrix = matrices[sector_index];

    // Extract subspace labels.
    int bra_L = bra_subspace.L();
    int bra_S = bra_subspace.S();
    int bra_J = bra_subspace.J();
    int bra_T = bra_subspace.T();
    int ket_L = ket_subspace.L();
    int ket_S = ket_subspace.S();
    int ket_J = ket_subspace.J();
    int ket_T = ket_subspace.T();

    if ((bra_T == ket_T) || (bra_S == ket_S)) {
      continue;
    }

    // Loop over bra and ket states.
#pragma omp parallel for collapse(2)
    for (std::size_t bra_index = 0; bra_index < bra_subspace.size();
         ++bra_index) {
      for (std::size_t ket_index = 0; ket_index < ket_subspace.size();
           ++ket_index) {
        const basis::RelativeStateLSJT bra_state(bra_subspace, bra_index);
        const basis::RelativeStateLSJT ket_state(ket_subspace, ket_index);

        // Extract state labels.
        int bra_n = bra_state.n();
        int ket_n = ket_state.n();

        // Common part of all radial integrals.
        Eigen::ArrayXd common_integrand = wt * scs_reg;
        common_integrand *=
            (ho_wfs.at(bra_L).row(bra_n) * ho_wfs.at(ket_L).row(ket_n));

        // Reduced matrix element calculation.
        double rme = 0;

        double tp_f =
            tp::CSpinTensorProductRME(bra_subspace, ket_subspace, 2, 1, 1);
        tp_f *= std::sqrt(10.);

        Eigen::ArrayXd y = common_integrand * zpir * ypir;
        y.head(1) = 0;  // Required to avoid divide by 0.
        y.tail(1) = 0;  // Required to avoid divide by 0.
        double integ_zpir_ypir = quadpp::spline::Integrate(x, y);

        rme = tp_f * integ_zpir_ypir;

        if (bra_L == ket_L) {
          double tp_g =
              tp::CSpinTensorProductRME(bra_subspace, ket_subspace, 0, 1, 1);

          y = common_integrand * tpir * ypir;
          y.head(1) = 0;  // Required to avoid divide by 0.
          y.tail(1) = 0;  // Required to avoid divide by 0.
          double integ_tpir_ypir = quadpp::spline::Integrate(x, y);

          rme += tp_g * integ_tpir_ypir;
        }

        rme *= tp::SpinTensorProductRME(bra_T, ket_T, 1);  // Isospin.
        rme *= -(mN * mPi * gA * gA) / (24 * constants::pi * FPi * FPi);

        matrix(bra_n, ket_n) = rme;
      }
    }
  }
}

}  // namespace relative
}  // namespace chime
