#ifndef IO_H
#define IO_H

#include <chrono>
#include <iomanip>
#include <array>
#include <boost/variant.hpp>
#include "basis/lsjt_operator.h"
#include "fmt/format.h"

namespace io
{
  using LSJTSpace = boost::variant<basis::RelativeSpaceLSJT,
                                   basis::RelativeCMSpaceLSJT>;
  using LSJTParameters = boost::variant<basis::RelativeOperatorParametersLSJT,
                                        basis::RelativeCMOperatorParametersLSJT>;
  using LSJTSectors = boost::variant<std::array<basis::RelativeSectorsLSJT, 3>,
                                     std::array<basis::RelativeCMSectorsLSJT, 3>>;

  inline void SetLSJTParameters(LSJTParameters& params,
                         basis::OperatorLabelsJT& op_labels,
                         int& Nmax, int& Jmax, bool& has_cm)
  {
    if (has_cm)
        params = basis::RelativeCMOperatorParametersLSJT(op_labels, Nmax);
    else
        params = basis::RelativeOperatorParametersLSJT(op_labels, Nmax, Jmax);
  }

  inline void SetSpace(LSJTSpace& space,
                          const basis::RelativeOperatorParametersLSJT& params)
  {
    space = basis::RelativeSpaceLSJT(params.Nmax, params.Jmax);
  }
  inline void SetSpace(LSJTSpace& space,
                          const basis::RelativeCMOperatorParametersLSJT& params)
  {
    space = basis::RelativeCMSpaceLSJT(params.Nmax);
  }

  template <class T, std::size_t size>
  std::array<T, size> GetArray()
  {
    std::array<T, size> a;
    return a;
  }

  inline void SetSectors(LSJTSectors& sectors, bool& has_cm)
  {
    if (has_cm)
      {
        sectors = GetArray<basis::RelativeCMSectorsLSJT, 3>();
      }
    else
      {
        sectors = GetArray<basis::RelativeSectorsLSJT, 3>();
      }
  }

  void PrintTruncationInfo(basis::RelativeOperatorParametersLSJT params)
  {
      fmt::print("Truncation: Nmax {0} Jmax {1} T0_max {2}\n",
                 params.Nmax, params.Jmax, params.T0_max);
  }

  void PrintTruncationInfo(basis::RelativeCMOperatorParametersLSJT params)
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
