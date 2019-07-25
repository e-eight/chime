#ifndef IO_H
#define IO_H

#include <chrono>
#include <iomanip>
#include <array>
#include <boost/variant.hpp>
#include "basis/lsjt_operator.h"
#include "fmt/format.h"
#include "constants.h"
#include "utility.h"

namespace io
{
  // Variant classes for holding both relative and relative-cm variants.
  using LSJTSpace = boost::variant<basis::RelativeSpaceLSJT,
                                   basis::RelativeCMSpaceLSJT>;
  using LSJTParameters = boost::variant<basis::RelativeOperatorParametersLSJT,
                                        basis::RelativeCMOperatorParametersLSJT>;
  using LSJTSectors = boost::variant<std::array<basis::RelativeSectorsLSJT, 3>,
                                     std::array<basis::RelativeCMSectorsLSJT, 3>>;
  using OscillatorLength = boost::variant<double, std::array<double, 2>>;

  // Choose the LSJTParameter type depending on whether cm is there or not.
  inline void SetLSJTParameters(LSJTParameters& params,
                         basis::OperatorLabelsJT& op_labels,
                         int& Nmax, int& Jmax, bool& has_cm)
  {
    if (has_cm)
        params = basis::RelativeCMOperatorParametersLSJT(op_labels, Nmax);
    else
        params = basis::RelativeOperatorParametersLSJT(op_labels, Nmax, Jmax);
  }

  // Choose the LSJTSpace type depending on whether cm is there or not.
  inline void SetLSJTSpace(LSJTSpace& space, int& Nmax, int& Jmax, bool& has_cm)
  {
    if (has_cm)
      space = basis::RelativeCMSpaceLSJT(Nmax, Jmax);
    else
      space = basis::RelativeSpaceLSJT(Nmax);
  }

  inline void SetLSJTSectors(LSJTSectors& sectors, bool& has_cm)
  {
    if (has_cm)
      {
        sectors = util::GetArray<basis::RelativeCMSectorsLSJT, 3>();
      }
    else
      {
        sectors = util::GetArray<basis::RelativeSectorsLSJT, 3>();
      }
  }

  inline void SetOscillatorLength(OscillatorLength& b, double& hw, bool& has_cm)
  {
    double b_rel = constants::hbarc
      / std::sqrt(constants::reduced_nucleon_mass_MeV * hw);
    if (!has_cm)
      {
        b = b_rel;
      }
    else
      {
        double b_cm = constants::hbarc
          / std::sqrt(2 * constants::nucleon_mass_MeV * hw);
        std::array<double, 2> b2 = {b_rel, b_cm};
        b = b2;
      }
  }

  void PrintTruncationInfo(basis::RelativeOperatorParametersLSJT& params)
  {
      fmt::print("Truncation: Nmax {0} Jmax {1} T0_max {2}\n",
                 params.Nmax, params.Jmax, params.T0_max);
  }

  void PrintTruncationInfo(basis::RelativeCMOperatorParametersLSJT& params)
  {
      fmt::print("Truncation: Nmax {0} Jmax {1} T0_max {2}\n",
                 params.Nmax, params.T0_max);
  }

  std::string GetTimeInfo()
  {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    auto time_str = std::to_string(time);
    return time_str;
  }

}

#endif
