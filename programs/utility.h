/*******************************************************************************
 util.h

 Defines auxiliary functions for use in chime.

 Language: C++11
 Soham Pal
 Iowa State University
*******************************************************************************/

#ifndef UTIL_H_
#define UTIL_H_

#include <cmath>
#include "constants.h"

namespace util {
  // Calculates the oscillator parameter length from the oscillator energy for the
  // relative coordinate.
  inline
  double brel(const double& oscillator_energy) {
    double result = constants::hbarc;
    result /= std::sqrt(constants::nucleon_mass_MeV * oscillator_energy / 2);
    return result;
  }

  // Calculates the oscillator parameter length from the oscillator energy for the
  // center of mass coordinate.
  inline
  double bcm(const double& oscillator_energy) {
    double result = constants::hbarc;
    result /= std::sqrt(2 * constants::nucleon_mass_MeV * oscillator_energy);
    return result;
  }
}

#endif
