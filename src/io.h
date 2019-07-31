#ifndef IO_H
#define IO_H

#include <chrono>
#include <iomanip>
#include <array>
#include <boost/variant.hpp>
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "fmt/format.h"
#include "chiral.h"
#include "constants.h"
#include "utility.h"

namespace io
{
  // Variant classes for holding both relative and relative-cm variants.
  // using LSJTSpace = boost::variant<basis::RelativeSpaceLSJT,
  //                                  basis::RelativeCMSpaceLSJT>;
  // using LSJTParameters = boost::variant<basis::RelativeOperatorParametersLSJT,
  //                                       basis::RelativeCMOperatorParametersLSJT>;
  // using LSJTSectors = boost::variant<std::array<basis::RelativeSectorsLSJT, 3>,
  //                                    std::array<basis::RelativeCMSectorsLSJT, 3>>;
  // using OscillatorLength = boost::variant<double, std::array<double, 2>>;

  // // Choose the LSJTParameter type depending on whether cm is there or not.
  // inline void SetLSJTParameters(LSJTParameters& params,
  //                        basis::OperatorLabelsJT& op_labels,
  //                        int& Nmax, int& Jmax, bool& has_cm)
  // {
  //   if (has_cm)
  //       params = basis::RelativeCMOperatorParametersLSJT(op_labels, Nmax);
  //   else
  //       params = basis::RelativeOperatorParametersLSJT(op_labels, Nmax, Jmax);
  // }

  // // Choose the LSJTSpace type depending on whether cm is there or not.
  // inline void SetLSJTSpace(LSJTSpace& space, int& Nmax, int& Jmax, bool& has_cm)
  // {
  //   if (has_cm)
  //     space = basis::RelativeCMSpaceLSJT(Nmax);
  //   else
  //     space = basis::RelativeSpaceLSJT(Nmax, Jmax);
  // }

  // inline void SetLSJTSectors(LSJTSectors& sectors, bool& has_cm)
  // {
  //   if (has_cm)
  //     {
  //       sectors = util::GetArray<basis::RelativeCMSectorsLSJT, 3>();
  //     }
  //   else
  //     {
  //       sectors = util::GetArray<basis::RelativeSectorsLSJT, 3>();
  //     }
  // }

  // inline void SetOscillatorLength(OscillatorLength& b, double& hw, bool& has_cm)
  // {
  //   double b_rel = constants::hbarc
  //     / std::sqrt(constants::reduced_nucleon_mass_MeV * hw);
  //   if (!has_cm)
  //     {
  //       b = b_rel;
  //     }
  //   else
  //     {
  //       double b_cm = constants::hbarc
  //         / std::sqrt(2 * constants::nucleon_mass_MeV * hw);
  //       std::array<double, 2> b2 = {b_rel, b_cm};
  //       b = b2;
  //     }
  // }

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

  using OperatorMatrix = std::array<basis::OperatorBlocks<double>, 3>;
  using RelativeSector = std::array<basis::RelativeSectorsLSJT, 3>;
  using RelativeCMSector = std::array<basis::RelativeCMSectorsLSJT, 3>;

  template <class StateType>
  void GenerateMatrixElements(std::size_t& T0,
                              std::size_t& sector_index,
                              std::size_t& bra_index,
                              std::size_t& ket_index,
                              StateType& bra,
                              StateType& ket,
                              OperatorMatrix& matrix,
                              OperatorMatrix& temp_matrix,
                              chiral::Operator& op,
                              chiral::Order& ord,
                              util::OscillatorParameter& b,
                              double& regulator,
                              bool& regularize)
  {
    auto rme = op.ReducedMatrixElement(ord, bra, ket, b, regulator, regularize);
    temp_matrix[T0][sector_index](bra_index, ket_index) = rme;
    matrix[T0][sector_index](bra_index, ket_index) += rme;
  }

  void GenerateRelativeMatrices(std::size_t& T0,
                                RelativeSector& sectors,
                                OperatorMatrix& matrix,
                                OperatorMatrix& temp_matrix,
                                chiral::Operator& op,
                                chiral::Order& ord,
                                util::OscillatorParameter& b,
                                double& regulator,
                                bool& regularize)
  {
    for (std::size_t sector_index = 0; sector_index < sectors[T0].size(); ++sector_index)
      {
        const basis::RelativeSectorsLSJT::SectorType& sector
          = sectors[T0].GetSector(sector_index);
        const basis::RelativeSectorsLSJT::SubspaceType& bra_subspace
          = sector.bra_subspace();
        const basis::RelativeSectorsLSJT::SubspaceType& ket_subspace
          = sector.ket_subspace();

        for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
          {
            const basis::RelativeStateLSJT bra_state(bra_subspace, bra_index);
            for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
              {
                const basis::RelativeStateLSJT ket_state(ket_subspace, ket_index);
                GenerateMatrixElements(T0, sector_index, bra_index, ket_index,
                                       bra_state, ket_state, matrix, temp_matrix,
                                       op, ord, b, regulator, regularize);
              }
          }
      }
  }

  void GenerateRelativeCMMatrices(std::size_t& T0,
                                  RelativeCMSector& sectors,
                                  OperatorMatrix& matrix,
                                  OperatorMatrix& temp_matrix,
                                  chiral::Operator& op,
                                  chiral::Order& ord,
                                  util::OscillatorParameter& b,
                                  double& regulator,
                                  bool& regularize)
  {
    for (std::size_t sector_index = 0; sector_index < sectors[T0].size(); ++sector_index)
      {
        const basis::RelativeCMSectorsLSJT::SectorType& sector
          = sectors[T0].GetSector(sector_index);
        const basis::RelativeCMSectorsLSJT::SubspaceType& bra_subspace
          = sector.bra_subspace();
        const basis::RelativeCMSectorsLSJT::SubspaceType& ket_subspace
          = sector.ket_subspace();

        for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
          {
            const basis::RelativeCMStateLSJT bra_state(bra_subspace, bra_index);
            for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
              {
                const basis::RelativeCMStateLSJT ket_state(ket_subspace, ket_index);
                GenerateMatrixElements(T0, sector_index, bra_index, ket_index,
                                       bra_state, ket_state, matrix, temp_matrix,
                                       op, ord, b, regulator, regularize);
              }
          }
      }
  }

  void WriteRelativeFiles(std::string& name,
                          std::string& order,
                          double& hw,
                          std::string& Nmax_str,
                          std::string& Jmax_str,
                          std::size_t& T0_min,
                          std::size_t& T0_max,
                          bool& regularize,
                          double& regulator,
                          basis::RelativeSpaceLSJT& space,
                          RelativeSector& sectors,
                          OperatorMatrix& matrix)
  {
    auto op = chiral::Operator::make(name);
    auto ord = chiral::Order::_from_string(order.c_str());
    util::OscillatorParameter b(hw);
    basis::OperatorLabelsJT labels(op->J0(), op->G0(), T0_min, T0_max,
                                   basis::SymmetryPhaseMode::kHermitian);
    basis::ConstructIdentityOperatorRelativeLSJT(labels, space, sectors, matrix);

    std::string hw_str = fmt::format("{:.1f}", hw);
    std::string time_str = GetTimeInfo();

    OperatorMatrix& temp_matrix = matrix;
    for (chiral::Order cord : chiral::Order::_values())
      {
        for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
          {
            GenerateRelativeMatrices(T0, sectors, matrix, temp_matrix,
                                     *op, ord, b, regulator, regularize);
          }

        std::string cord_str = cord._to_string();
        std::string order_filename
          = (name + "_2b_rel_" + cord_str + "_N" + Nmax_str + "_J" + Jmax_str
             + "_hw" + hw_str + "_" + time_str + ".dat");
        basis::WriteRelativeOperatorLSJT(order_filename, space, labels,
                                         sectors, temp_matrix, true);

        if (cord == ord)
          break;
      }
    std::string ord_str = ord._to_string();
    std::string cumulative_filename
      = (name + "_2b_rel_" + ord_str + "_cumulative" + "_N" + Nmax_str
         + "_J" + Jmax_str + hw_str + "_" + time_str + ".dat");
    basis::WriteRelativeOperatorLSJT(cumulative_filename, space,
                                     labels, sectors, matrix, true);
  }

  void WriteRelativeCMFiles(std::string& name,
                            std::string& order,
                            double& hw,
                            std::string& Nmax_str,
                            std::size_t& T0_min,
                            std::size_t& T0_max,
                            bool& regularize,
                            double& regulator,
                            basis::RelativeCMSpaceLSJT& space,
                            RelativeCMSector& sectors,
                            OperatorMatrix& matrix)
  {
    auto op = chiral::Operator::make(name);
    auto ord = chiral::Order::_from_string(order.c_str());
    util::OscillatorParameter b(hw);
    basis::OperatorLabelsJT labels(op->J0(), op->G0(), T0_min, T0_max,
                                   basis::SymmetryPhaseMode::kHermitian);

    std::string hw_str = fmt::format("{:.1f}", hw);
    std::string time_str = GetTimeInfo();

    OperatorMatrix& temp_matrix = matrix;
    for (chiral::Order cord : chiral::Order::_values())
      {
        for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
          {
            GenerateRelativeCMMatrices(T0, sectors, matrix, temp_matrix,
                                       *op, ord, b, regulator, regularize);
          }

        std::string cord_str = cord._to_string();
        std::string order_filename
          = (name + "_2b_rel_" + cord_str + "_N" + Nmax_str
             + "_hw" + hw_str + "_" + time_str + ".dat");
        basis::WriteRelativeCMOperatorLSJT(order_filename, space, labels,
                                           sectors, temp_matrix, true);

        if (cord == ord)
          break;
      }
    std::string ord_str = ord._to_string();
    std::string cumulative_filename
      = (name + "_2b_rel_" + ord_str + "_cumulative" + "_N" + Nmax_str
         + hw_str + "_" + time_str + ".dat");
    basis::WriteRelativeCMOperatorLSJT(cumulative_filename, space,
                                     labels, sectors, matrix, true);
  }

  void WriteFiles(std::string& name,
                  std::string& order,
                  bool& has_cm,
                  double& hw,
                  std::size_t& Nmax,
                  std::size_t& Jmax,
                  std::size_t& T0_min,
                  std::size_t& T0_max,
                  bool& regularize,
                  double& regulator)
  {
    std::string Nmax_str = fmt::format("{:d}", Nmax);
    std::string Jmax_str = fmt::format("{:d}", Jmax);
    if (!has_cm)
      {
        basis::RelativeSpaceLSJT space(Nmax, Jmax);
        RelativeSector sectors;
        OperatorMatrix matrix;
        WriteRelativeFiles(name, order, hw, Nmax_str, Jmax_str, T0_min, T0_max,
                           regularize, regulator, space, sectors, matrix);
      }
    else
      {
        basis::RelativeCMSpaceLSJT space(Nmax);
        RelativeCMSector sectors;
        OperatorMatrix matrix;
        WriteRelativeCMFiles(name, order, hw, Nmax_str, T0_min, T0_max,
                             regularize, regulator, space, sectors, matrix);
      }
  }

}

#endif
