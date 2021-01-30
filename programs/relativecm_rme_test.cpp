#include "relativecm_rme.h"

#include <fstream>

#include "basis_func/ho.h"
#include "chime.h"
#include "constants.h"
#include "quadpp/quadpp.h"
#include "quadpp/spline.h"

int main()
{
  int Nmax = 20;
  double hbomega = 20;
  double R = 1.;
  double mPi = chime::constants::pion_mass_fm;

  int npts = 1001;
  Eigen::ArrayXd x, r, jac;
  quadpp::SemiInfiniteIntegralMesh(npts, 0, 1, x, r, jac);
  Eigen::ArrayXd wt = r * r * jac;

  Eigen::ArrayXXd psi(Nmax + 1, npts);
  double brel = chime::RelativeOscillatorLength(hbomega);
  double bcm = chime::CMOscillatorLength(hbomega);

  basis_func::ho::WaveFunctionsUptoMaxN(psi, r, Nmax, Nmax, brel,
                                        basis_func::Space::coordinate);
  // psi.col(npts - 1) = 0;
  // std::cout << psi << "\n";

  // Relative radial integral kernels.
  Eigen::ArrayXd mpir = mPi * r;
  // mpir.tail(1) = 0;
  Eigen::ArrayXd expmpir = Eigen::exp(-mpir);
  // expmpir.tail(1) = 0;
  // expmpir.head(1) = 0;
  Eigen::ArrayXd ypir = expmpir / mpir;
  // ypir.head(1) = 0;
  // ypir.tail(1) = 0;// This needs to be set to avoid divide by 0;
  Eigen::ArrayXd zpir = (1. + mpir);
  // zpir.head(1) = 0;
  // zpir.tail(1) = 0;
  Eigen::ArrayXd tpir = (-1. + 2 * mpir);
  // tpir.head(1) = 0;
  // tpir.tail(1) = 0;
  Eigen::ArrayXd wpir = (1. + (3 * zpir / mpir.square()));
  // wpir.head(1) = 0; // This needs to be set to avoid divide by 0;
  // wpir.tail(1) = 0;

  // Semilocal coordinate space regulator.
  Eigen::ArrayXd scs_reg = Eigen::ArrayXd::Ones(npts);
  if (R != 0) {
    scs_reg *= Eigen::pow(1. - Eigen::exp(-(r * r) / (R * R)), 6);
  }
  // scs_reg.head(1) = 0;
  // scs_reg.tail(1) = 1;

  // std::cout << "mpir" << "\n";
  // std::cout << mpir << "\n";
  // std::cout << "expmpir" << "\n";
  // std::cout << expmpir << "\n";
  // std::cout << "ypir" << "\n";
  // std::cout << ypir << "\n";
  // std::cout << "zpir" << "\n";
  // std::cout << zpir << "\n";
  // std::cout << "tpir" << "\n";
  // std::cout << tpir << "\n";
  // std::cout << "wpir" << "\n";
  // std::cout << wpir << "\n";
  // std::cout << "scs_reg" << "\n";
  // std::cout << scs_reg << "\n";

  for (int bra_nr = 0; bra_nr <= Nmax; ++bra_nr) {
    for (int ket_nr = 0; ket_nr <= bra_nr; ++ket_nr) {
      Eigen::ArrayXd psin2 = psi.row(bra_nr) * psi.row(ket_nr);
      Eigen::ArrayXd y = wt * scs_reg * psin2 * zpir * ypir;

      y.head(1) = 0;
      y.tail(1) = 0;
      double integ_zpir_ypir = quadpp::spline::Integrate(x, y);

      std::cout << bra_nr << " " << ket_nr << " " << integ_zpir_ypir << "\n";
    }
  }

  // Eigen::ArrayXXd output(npts-1, 2);
  // output.col(0) = r.head(npts-1);
  // Eigen::ArrayXd psisq = psi.row(Nmax) * psi.row(Nmax);
  // output.col(1) = y.head(npts-1);
  // Eigen::IOFormat fmt(7, Eigen::DontAlignCols, " ", "\n");

  // std::ofstream file;
  // file.open("integrand.txt");
  // file << output.format(fmt) << "\n";
  // file.close();
}
