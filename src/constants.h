#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_math.h>

namespace constants
{
  constexpr double pi = M_PI;
  constexpr double sqrtpi = M_SQRTPI;
  constexpr double sqrt2 = M_SQRT2;

  // Conversion Factor
  constexpr double hbarc = 197.3269602; // (in MeV fm)
  constexpr double hbarc_GeV = 0.1973269602; // (in GeV fm)

  // Masses
  constexpr double pion_mass_MeV = 134.9770;
  constexpr double pion_mass_fm = pion_mass_MeV / hbarc; // (in fm^{-1})
  constexpr double reduced_nucleon_mass_MeV = 469.4593340;
  constexpr double reduced_nucleon_mass_fm = reduced_nucleon_mass_MeV / hbarc;
  constexpr double nucleon_mass_MeV = reduced_nucleon_mass_MeV * 2;
  constexpr double nucleon_mass_fm = reduced_nucleon_mass_fm * 2;
  constexpr double nuclear_magneton_MeV = 1.0 / (nucleon_mass_MeV * 2); // (e = 1)
  constexpr double nuclear_magneton_fm = 1.0 / (nucleon_mass_fm * 2); // (e = 1)

  // Low energy constants
  constexpr double isoscalar_nucleon_charge_radius_sq_fm = 0.603729; // (in fm^2)
  constexpr double isoscalar_nucleon_magnetic_moment = 0.88;
  constexpr double gA = 1.25;
  constexpr double pion_decay_constant_MeV = 92.4;
  constexpr double pion_decay_constant_fm = pion_decay_constant_MeV / hbarc; // (in fm^{-1})
  constexpr double d9_GeV = -0.011; // (in GeV^{-2})
  constexpr double d9_fm = d9_GeV * hbarc_GeV * hbarc_GeV; // (in fm^2)
  constexpr double L2_GeV = 0.188; // (in GeV^{-4})
  constexpr double L2_fm = L2_GeV * hbarc_GeV * hbarc_GeV * hbarc_GeV * hbarc_GeV; // (in fm^4)
}

#endif
