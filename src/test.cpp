#include "rme_extras.h"
#include "integrals.h"
#include "threedho.h"
#include "utility.h"
#include "constants.h"
#include "fmt/format.h"
#include "chiral.h"

int main()
{
  fmt::print("Spin Symmetric: \n");
  for (int i = 0; i <= 1; ++i)
    {
      for (int j = 0; j <=1; ++j)
        {
          fmt::print("{} {} {}\n", i, j, am::SpinSymmetricRME(i, j));
        }
    }

  fmt::print("Spin Antisymmetric: \n");
  for (int i = 0; i <= 1; ++i)
    {
      for (int j = 0; j <=1; ++j)
        {
          fmt::print(" {}\n", am::SpinAntisymmetricRME(i, j));
        }
    }

  fmt::print("Relative Spin Antisymmetric: \n");
  for (int l = 0; l <= 2; l++)
    {
      fmt::print(" {}\n", am::RelativeSpinAntisymmetricRME(l, l, 1, 0, 1, 0, 0, 1));
    }

  fmt::print("LSCoupledCRME: \n");
  for (int l = 0; l <= 2; l++)
    {
      fmt::print(" {}\n", am::SphericalHarmonicCRME(l, l, 0));
    }

  fmt::print("Wigner 6J: {}\n", am::Wigner6J(0, 0, 0, 1, 1, 1));

  fmt::print("Wigner 9J: {}\n", am::Wigner9J(1, 1, 1, 0, 1, 1, 1, 1, 0));

  fmt::print("Relative spin symmetric product: {}\n",
             am::RelativeSpinSymmetricRME(1, 1, 1, 1, 0, 1, 0, 1));

  fmt::print("Relative spin symmetric direct: {}\n", am::Wigner9J(1, 1, 1, 0, 1, 1, 1, 1, 0) * HatProduct(1, 1, 1, 1) * am::SphericalHarmonicCRME(1, 1, 0) * am::SpinSymmetricRME(1, 1));

  fmt::print("C rme: {}\n", am::SphericalHarmonicCRME(1, 1, 0));

  fmt::print("Relative spin isospin asymmetric product: {}\n",
             am::RelativeSpinAntisymmetricRME(0, 0, 1, 0, 1, 0, 0, 1)
             * am::SpinAntisymmetricRME(0, 1));

  fmt::print("Relative Lrel: {} {} {}\n",
             am::RelativeLrelRME(1, 1, 1, 1, 0, 1),
             am::RelativeLrelRME(1, 1, 1, 1, 0, 1),
             am::RelativeLrelRME(2, 2, 1, 1, 3, 3));

  fmt::print("LECs: D = {}\n", constants::D_fm);

  fmt::print("Hat product: {}\n", HatProduct(1, 1, 1, 1));

  fmt::print("(Delta regularization) F: n=0: {}, n=100: {}\n",
             quadrature::F(0, 0, 1), quadrature::F(100, 0, 1));

  fmt::print("Chiral orders: ");
  for(const auto& order : chiral::v_order)
    fmt::print("{} ", chiral::reverse_m_order[order]);
  fmt::print("\n");

  auto num = constants::gA * cube(constants::pion_mass_fm) * constants::d18_fm;
  auto denom = 12 * constants::pi * square(constants::pion_decay_constant_fm);
  auto lec_prefactor = num / denom;
  lec_prefactor /= constants::nuclear_magneton_fm;
  fmt::print("d18 product: {}\n", lec_prefactor);
}
