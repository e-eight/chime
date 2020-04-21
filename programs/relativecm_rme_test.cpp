#include <fstream>
#include "constants.h"
#include "relativecm_rme.h"

int main()
{
  int Nmax = 20;
  double hbomega = 20;
  double R = 1.;
  double mPi = constants::pion_mass_fm;

  int npts = 1001;
  Eigen::ArrayXd x = Eigen::ArrayXd::LinSpaced(npts, 0., 1.);
  Eigen::ArrayXd r = x / (1. - x);
  Eigen::ArrayXd dr = 1. / Eigen::square(1. - x);
  Eigen::ArrayXd wt = r * r * dr;

  Eigen::ArrayXXd psi(Nmax + 1, npts);
  double brel = util::brel(hbomega);
  double bcm = util::bcm(hbomega);

  basis_func::ho::WF(psi, r, Nmax, Nmax, brel, "coordinate");
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
      double integ_zpir_ypir = quad::spline::Integrate(x, y);

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
