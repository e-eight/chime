/*******************************************************************************
 chime.h

 Utility functions and classes for chime.

 Language: C++11
 Soham Pal
 Iowa State University
*******************************************************************************/

#ifndef CHIME_H_
#define CHIME_H_

#include <cmath>
#include <fstream>
#include "basis/jt_operator.h"
#include "constants.h"

namespace chime {
// Calculates the oscillator parameter length from the oscillator energy for
// the relative coordinate.
inline double RelativeOscillatorLength(const double& oscillator_energy)
{
  double result = constants::hbarc;
  result /= std::sqrt(constants::nucleon_mass_MeV * oscillator_energy / 2);
  return result;
}

// Calculates the oscillator parameter length from the oscillator energy for
// the center of mass coordinate.
inline double CMOscillatorLength(const double& oscillator_energy)
{
  double result = constants::hbarc;
  result /= std::sqrt(2 * constants::nucleon_mass_MeV * oscillator_energy);
  return result;
}

// Calculates the LENPIC Semilocal Coordinate Space (SCS) regulator.
inline double SCSRegulator(const double& r, const double& R)
{
  if (R == 0) {
    return 1.0;
  }
  double result = std::pow(1 - std::exp(-(r * r) / (R * R)), 6);
  return result;
}

}  // namespace chime

#endif
