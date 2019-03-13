#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <gsl/gsl_math.h>

namespace constants
{
  constexpr double pi = M_PI;

  // Conversion Factor
  constexpr double hbarc = 197.3269602; // (in MeV fm)

  // Masses
  constexpr double pion_mass_MeV = 134.9770;
  constexpr double pion_mass_fm = pion_mass_MeV / hbarc; // (in fm^{-1})
  constexpr double reduced_nucleon_mass_MeV = 469.4593340;
  constexpr double reduced_nucleon_mass_fm = reduced_nucleon_mass_MeV / hbarc;
  constexpr double nucleon_mass_MeV = reduced_nucleon_mass_MeV * 2;
  constexpr double nucleon_mass_fm = reduced_nucleon_mass_fm * 2;

  // Low energy constants
  constexpr double isoscalar_nucleon_charge_radius_sq_fm = 0.603729; // (in fm^2)
  constexpr double gA = 1.25;
  constexpr double pion_decay_constant_MeV = 92.4;
  constexpr double pion_decay_constant_fm = pion_decay_constant_MeV / hbarc; // (in fm^{-1})
}

#endif
