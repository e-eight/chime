#include "relative_rme.h"

#include <fstream>

#include "basis_func/ho.h"
#include "chime.h"
#include "constants.h"
#include "quadpp/quadpp.h"
#include "quadpp/spline.h"
#include "tprme.h"

int main()
{
  constexpr double mPi = chime::constants::pion_mass_fm;
  constexpr double mN = chime::constants::nucleon_mass_fm;
  constexpr double FPi = chime::constants::pion_decay_constant_fm;
  constexpr double gA = chime::constants::gA;

  int Lp = 0, Sp = 0, Jp = 0, Tp = 1;
  int L = 0, S = 1, J = 1, T = 0;
  int Nmax = 2, g = 0;
  int Np = 16, N = 16;
  int np = (Np - Lp) / 2, n = (N - L) / 2;
  basis::RelativeSubspaceLSJT bra_subspace(Lp, Sp, Jp, Tp, g, Nmax);
  basis::RelativeSubspaceLSJT ket_subspace(L, S, J, T, g, Nmax);

  double tp_f =
      chime::tp::CSpinTensorProductRME(bra_subspace, ket_subspace, 2, 1, 1);
  double tp_g =
      chime::tp::CSpinTensorProductRME(bra_subspace, ket_subspace, 0, 1, 1);
  std::cout << "tp_f: " << tp_f << " tp_g: " << tp_g << "\n";

  std::size_t npts = 1001;
  Eigen::ArrayXd x, r, jac;
  quadpp::SemiInfiniteIntegralMesh(npts, 0, 1, x, r, jac);
  Eigen::ArrayXd wt = r * r * jac;

  Eigen::ArrayXd mpir = mPi * r;
  Eigen::ArrayXd expmpir = Eigen::exp(-mpir);
  Eigen::ArrayXd ypir = expmpir / mpir;
  Eigen::ArrayXd zpir = (1. + mpir);
  Eigen::ArrayXd tpir = (-1. + 2 * mpir);

  Eigen::ArrayXd scs_reg = Eigen::pow(1. - Eigen::exp(-(r * r)), 6);

  double b = chime::RelativeOscillatorLength(20);
  Eigen::ArrayXXd bra_wfs(Nmax + 1, r.size()), ket_wfs(Nmax + 1, r.size());
  basis_func::ho::WaveFunctionsUptoMaxN(bra_wfs, r, Nmax, Lp, b,
                                        basis_func::Space::coordinate);
  basis_func::ho::WaveFunctionsUptoMaxN(ket_wfs, r, Nmax, L, b,
                                        basis_func::Space::coordinate);

  Eigen::ArrayXd y;
  y = tpir * ypir * scs_reg * wt;
  y *= bra_wfs.row(np) * ket_wfs.row(n);
  y.head(1) = 0;
  y.tail(1) = 0;
  double integ_tpi_ypi = quadpp::spline::Integrate(x, y);
  std::cout << "integ_tpi_ypi: " << integ_tpi_ypi << " \n";

  std::cout << "tp_g * integ_tpi_ypi: " << tp_g * integ_tpi_ypi << "\n";

  double prefactor =
      -(mPi * mN * gA * gA) / (24 * chime::constants::pi * FPi * FPi);
  double isospin = chime::tp::SpinTensorProductRME(Tp, T, 1);
  std::cout << "prefactor: " << prefactor << "\n";
  std::cout << "isospin: " << isospin << "\n";
  std::cout << "tp_g full result: "
            << prefactor * isospin * tp_g * integ_tpi_ypi << "\n";
}
