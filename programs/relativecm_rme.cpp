#include "relativecm_rme.h"

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
namespace relcm {

constexpr double mPi = constants::pion_mass_fm;
constexpr double mN = constants::nucleon_mass_fm;
constexpr double FPi = constants::pion_decay_constant_fm;
constexpr double gA = constants::gA;

void ConstructMu2nNLOOperator(
    const basis::RelativeCMOperatorParametersLSJT& op_params,
    const basis::RelativeCMSpaceLSJT& relcm_space,
    std::array<basis::RelativeCMSectorsLSJT, 3>& relcm_sectors,
    std::array<basis::OperatorBlocks<double>, 3>& relcm_matrices,
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

  const int npts = 3001;
  const double low = 0, high = 1;
  Eigen::ArrayXd x, r, jac, wt;
  quadpp::SemiInfiniteIntegralMesh(npts, low, high, x, r, jac);
  wt = r * r * jac;  // weights for radial integral with transformed variable

  std::cout << "  Generating basis functions...\n";
  std::vector<Eigen::ArrayXXd> ho_wfs;
  double brel = chime::RelativeOscillatorLength(oscillator_energy);
  double bcm = chime::CMOscillatorLength(oscillator_energy);
  basis_func::ho::WaveFunctionsUptoMaxL(ho_wfs, r, Nmax, Nmax, brel,
                                        basis_func::Space::coordinate);

  // Relative radial integral kernels.
  std::cout << "  Generating integration kernels...\n";
  Eigen::ArrayXd mpir = mPi * r;
  Eigen::ArrayXd expmpir = Eigen::exp(-mpir);
  Eigen::ArrayXd ypir = expmpir / mpir;
  Eigen::ArrayXd zpir = (1. + mpir);
  Eigen::ArrayXd tpir = (-1. + 2 * mpir);
  Eigen::ArrayXd wpir = (1. + (3 * zpir / mpir.square()));

  // Semilocal coordinate space regulator.
  Eigen::ArrayXd scs_reg = Eigen::ArrayXd::Ones(npts);
  if (R != 0) {
    scs_reg *= Eigen::pow(1. - Eigen::exp(-(r * r) / (R * R)), 6);
  }

  // Zero initialize operator.
  std::cout << "  Zero initializing operator...\n";
  for (int T = op_params.T0_min; T <= op_params.T0_max; ++T) {
    relcm_sectors[T] = basis::RelativeCMSectorsLSJT(relcm_space, op_params.J0,
                                                    T, op_params.g0);
    basis::SetOperatorToZero(relcm_sectors[T], relcm_matrices[T]);
  }

  // Select T0 component.
  const basis::RelativeCMSectorsLSJT& sectors = relcm_sectors[T0];
  basis::OperatorBlocks<double>& matrices = relcm_matrices[T0];

  // Reduced matrix element calculation.
  std::cout << "  Starting matrix element calculation...\n";
  for (std::size_t sector_index = 0; sector_index < sectors.size();
       ++sector_index) {
    const basis::RelativeCMSectorsLSJT::SectorType& sector =
        sectors.GetSector(sector_index);
    const basis::RelativeCMSubspaceLSJT& bra_subspace = sector.bra_subspace();
    const basis::RelativeCMSubspaceLSJT& ket_subspace = sector.ket_subspace();

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

    if (bra_T == ket_T) {
      continue;
    }

    // Loop over bra and ket states.
#pragma omp parallel for collapse(2)
    for (std::size_t bra_index = 0; bra_index < bra_subspace.size();
         ++bra_index) {
      for (std::size_t ket_index = 0; ket_index < ket_subspace.size();
           ++ket_index) {
        const basis::RelativeCMStateLSJT bra_state(bra_subspace, bra_index);
        const basis::RelativeCMStateLSJT ket_state(ket_subspace, ket_index);

        // Extract state labels.
        int bra_nr = bra_state.Nr();
        int bra_lr = bra_state.lr();
        int bra_nc = bra_state.Nc();
        int bra_lc = bra_state.lc();
        int ket_nr = ket_state.Nr();
        int ket_lr = ket_state.lr();
        int ket_nc = ket_state.Nc();
        int ket_lc = ket_state.lc();

        // Common part of all radial integrals.
        Eigen::ArrayXd common_integrand = wt * scs_reg;
        common_integrand *=
            ho_wfs.at(bra_lr).row(bra_nr) * ho_wfs.at(ket_lr).row(ket_nr);

        // Reduced matrix element calculation.
        // Pauli matrix tensor product in spin space enforces the bra and
        // ket spins to be the same for the relative-cm part, and to be
        // different for the purely relative part.
        double rme = 0;

        if (bra_S == ket_S) {
          // Relative-cm part.
          double tp_a =
              tp::CCSpinTensorProductRME(bra_state, ket_state, 1, 1, 1, 0, 1);
          tp_a *= -std::sqrt(3.);

          Eigen::ArrayXd y = common_integrand * expmpir;
          y.head(1) = 0;  // Required to avoid divide by 0.
          y.tail(1) = 0;  // Required to avoid divide by 0.
          double integ_expmpir = quadpp::spline::Integrate(x, y);

          rme = tp_a * integ_expmpir;

          if (bra_S == 1) {
            // Rank 2 Pauli Matrix tensor product.

            double tp_b =
                tp::CCSpinTensorProductRME(bra_state, ket_state, 1, 1, 1, 2, 1);
            tp_b *= std::sqrt(3. / 5.);

            double tp_c =
                tp::CCSpinTensorProductRME(bra_state, ket_state, 1, 1, 2, 2, 1);
            tp_c *= std::sqrt(9. / 5.);

            double tp_d =
                tp::CCSpinTensorProductRME(bra_state, ket_state, 3, 1, 2, 2, 1);
            tp_d *= std::sqrt(14. / 5.);

            double tp_e =
                tp::CCSpinTensorProductRME(bra_state, ket_state, 3, 1, 3, 2, 1);
            tp_e *= std::sqrt(28. / 5.);

            y *= wpir;
            y.head(1) = 0;  // Required to avoid divide by 0.
            y.tail(1) = 0;  // Required to avoid divide by 0.
            double integ_expmpir_wpir = quadpp::spline::Integrate(x, y);

            rme += (tp_b + tp_c + tp_d + tp_e) * integ_expmpir_wpir;
          }

          double integ_cm = 0;  // CM coordinate integral; analytical result.
          if (bra_lc == ket_lc + 1) {
            integ_cm = ((std::sqrt(ket_nc + ket_lc + 1.5) * (bra_nc == ket_nc))
                        + (std::sqrt(ket_nc) * (bra_nc + 1 == ket_nc)));
          }
          else if (bra_lc + 1 == ket_lc) {
            integ_cm = ((std::sqrt(bra_nc + ket_nc + 1.5)) * (bra_nc == ket_nc)
                        + (std::sqrt(bra_nc) * (bra_nc == ket_nc + 1)));
          }
          integ_cm *= mPi * bcm;

          rme *= integ_cm;
        }
        else {
          // Purely relative part. The cm labels for the bra and ket
          // must be the same.
          if ((bra_nc == ket_nc) && (bra_lc == ket_lc)) {
            double tp_f =
                tp::CCSpinTensorProductRME(bra_state, ket_state, 2, 0, 2, 1, 1);
            tp_f *= std::sqrt(10.);

            Eigen::ArrayXd y = common_integrand * zpir * ypir;
            y.head(1) = 0;  // Required to avoid divide by 0.
            y.tail(1) = 0;  // Required to avoid divide by 0.
            double integ_zpir_ypir = quadpp::spline::Integrate(x, y);

            rme = tp_f * integ_zpir_ypir;

            if (bra_lr == ket_lr) {
              // Rank 0 spherical harmonic.
              double tp_g = tp::CCSpinTensorProductRME(bra_state, ket_state, 0,
                                                       0, 0, 1, 1);

              y = common_integrand * tpir * ypir;
              y.head(1) = 0;  // Required to avoid divide by 0.
              y.tail(1) = 0;  // Required to avoid divide by 0.
              double integ_tpir_ypir = quadpp::spline::Integrate(x, y);

              rme += tp_g * integ_tpir_ypir;
            }
          }
        }
        rme *= tp::SpinTensorProductRME(bra_T, ket_T, 1);  // Isospin.
        rme *= -(mN * mPi * gA * gA) / (24 * constants::pi * FPi * FPi);

        matrix(bra_index, ket_index) = rme;
      }
    }
  }
}

}  // namespace relcm
}  // namespace chime
