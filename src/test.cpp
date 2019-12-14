#include "rme_extras.h"
#include "integrals.h"
#include "threedho.h"
#include "utility.h"
#include "constants.h"
#include "fmt/format.h"
#include "chiral.h"
#include "m1.h"
// #include "tprme.h"

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

  fmt::print("Chiral orders: ");
  for(const auto& order : chiral::v_order)
    fmt::print("{} ", chiral::reverse_m_order[order]);
  fmt::print("\n");

  auto nmass = constants::nucleon_mass_MeV;
  auto rnmass = constants::reduced_nucleon_mass_MeV;
  fmt::print("Average nucleon mass: {}\n", nmass);
  fmt::print("Reduced nucleon mass: {}\n", rnmass);
  fmt::print("mass ratio: {}\n", nmass / rnmass);

  auto cnorm = ho::CoordinateSpaceNorm(2, 2, 1);
  fmt::print("Coordinate space HO wf norm (n=2, l=2, b=1): {}\n", cnorm);

  fmt::print("Yukawa integrals: \n");
  auto b = ho::OscillatorParameter(20);
  auto brel = b.relative();
  auto bcm = b.cm();
  auto mpi = constants::pion_mass_fm;
  fmt::print("ZπYπ: \n");
  fmt::print("{}\n", quadrature::IntegralZPiYPiR(0, 0, 0, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralZPiYPiR(0, 0, 1, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralZPiYPiR(1, 1, 1, 1, brel, mpi, true, 1.0));
  fmt::print("TπYπ: \n");
  fmt::print("{}\n", quadrature::IntegralTPiYPiR(0, 0, 0, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralTPiYPiR(0, 0, 1, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralTPiYPiR(1, 0, 1, 0, brel, mpi, true, 1.0));
  fmt::print("(mπ r)Yπ: \n");
  fmt::print("{}\n", quadrature::IntegralMPiRYPiR(0, 0, 0, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralMPiRYPiR(0, 1, 0, 0, brel, mpi, true, 1.0));
  fmt::print("(mπ r)WπYπ: \n");
  fmt::print("{}\n", quadrature::IntegralMPiRWPiRYPiR(0, 0, 0, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralMPiRWPiRYPiR(0, 1, 0, 0, brel, mpi, true, 1.0));
  fmt::print("{}\n", quadrature::IntegralMPiRWPiRYPiR(0, 3, 0, 0, brel, mpi, true, 1.0));

  auto gpi = constants::gA / (2 * constants::pion_decay_constant_fm);
  auto mN = constants::nucleon_mass_fm;
  auto prefactor = -mpi * mN * square(gpi) / (3 * constants::pi);
  fmt::print("M1 two-body prefactor: {}\n", prefactor);

  fmt::print("U1...\n");
  fmt::print("a: {} {} {} \n",
             chiral::U1aRME(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             chiral::U1aRME(1, 1, 1, 0, 1, 0, 0, 0, 0, 0),
             chiral::U1aRME(0, 0, 0, 1, 1, 1, 1, 1, 1, 1));
  fmt::print("bcde: {} {} {} \n",
             chiral::U1bcdeRME(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
             chiral::U1bcdeRME(1, 1, 1, 1, 1, 0, 0, 0, 1, 1),
             chiral::U1bcdeRME(0, 0, 0, 1, 1, 1, 1, 1, 1, 1));
  fmt::print("f: {} {} {} \n",
             chiral::U1fS1RME(0, 0, 0, 0, 0, 1, 0, 1, 1, 1),
             chiral::U1fS1RME(0, 0, 0, 1, 1, 1, 0, 1, 0, 1),
             chiral::U1fS1RME(1, 1, 0, 0, 0, 1, 1, 2, 1, 1));
  fmt::print("S1: {} \n", chiral::S1RME(1, 1, 0, 0, 0, 1, 1, 0, 1, 1));

  // double u1a_tprme = (-std::sqrt(3)
  //                     * am::tp::RelativeCMSpinRME(0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, am::tp::SpinProductRME, 0, 1));
  // fmt::print("LS Coupled TPRME: {} \n", u1a_tprme);

}
