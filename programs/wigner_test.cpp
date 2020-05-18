#include "am/am.h"
#include "am/wigner_gsl.h"
#include "tprme.h"

int main()
{
  std::cout << am::Wigner9J(0, 1, 1, 0, 1, 1, 0, 0, 0) << "\n";

  int lp = 0, sp = 0, jp = 0, tp = 1;
  int l = 2, s = 1, j = 1, t = 0;
  int Nmax = 2, g = 0;

  std::cout << "C2S1 with Wigner 9J formula \n";
  basis::RelativeSubspaceLSJT bra_subspace(lp, sp, jp, tp, g, Nmax);
  basis::RelativeSubspaceLSJT ket_subspace(l, s, j, t, g, Nmax);
  double C2S1n =
      chime::tp::CSpinTensorProductRME(bra_subspace, ket_subspace, 2, 1, 1);
  std::cout << C2S1n << "\n";

  std::cout << "C2S1 with Wigner 6J formula \n";
  double C2S1s = ParitySign(lp + l) * am::Wigner6J(2, 1, 1, j, l, lp);
  C2S1s *= Hat(j) * std::sqrt(6) * am::SphericalHarmonicCRME(lp, l, 2);
  std::cout << C2S1s << "\n";
}
