#ifndef IO_H
#define IO_H

#include <chrono>
#include <iomanip>
#include <array>
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "fmt/format.h"
#include "chiral.h"
#include "utility.h"
#include "operators.h"

namespace io
{
  std::string GetTimeInfo()
  {
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);
    auto time_str = std::to_string(time);
    return time_str;
  }

  void WriteRelativeFiles(const std::string& name,
                          const std::string& order,
                          const std::size_t& Abody,
                          const double& hw,
                          const std::size_t& Nmax,
                          const std::size_t& Jmax,
                          const std::size_t& T0_min,
                          const std::size_t& T0_max,
                          const bool& regularize,
                          const double& regulator)
  {
    // Create operator.
    fmt::print("Generating {} matrix elements...\n", name);
    auto op = chiral::Operator::make(name);
    chiral::Order ord;
    // try
    //   {
    //     ord = chiral::m_order.at(order);
    //   }
    // catch(std::out_of_range& oor)
    //   {
    //     std::cerr << "Chiral order must be in [lo, n4lo].\n";
    //     std::terminate();
    //   }
    ord = chiral::m_order.at(order);

    // Set up operator and file header parameters.
    fmt::print("Beginning RelativeLSJT operator basis setup...\n");
    basis::OperatorLabelsJT labels(op->J0(), op->G0(), T0_min, T0_max,
                                   basis::SymmetryPhaseMode::kHermitian);
    basis::RelativeOperatorParametersLSJT params(labels, Nmax, Jmax);

    // Set up relative space.
    basis::RelativeSpaceLSJT space(Nmax, Jmax);

    // Set up operator containers.
    // These arrays store the information for each isospin channel.
    std::array<basis::RelativeSectorsLSJT, 3> sectors;
    std::array<basis::OperatorBlocks<double>, 3> matrices;

    // Populate containers.
    basis::ConstructZeroOperatorRelativeLSJT(labels, space, sectors, matrices);

    // Print truncation information.
    fmt::print("Truncation: Nmax {:d} Jmax {:d} T0_max {:d}\n", Nmax, Jmax, T0_max);
    for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
      {
        fmt::print(" T0 {:d} Upper triangular entries {:d} Allocated entries {:d}\n",
                   T0, basis::UpperTriangularEntries(sectors[T0]),
                   basis::AllocatedEntries(matrices[T0]));
      }

    // Time information.
    std::string time_str = GetTimeInfo();

    // Temporary matrices for storing contributions at each order.
    std::array<basis::OperatorBlocks<double>, 3> temp_matrices = matrices;

    // Iterate over chiral orders.
    for (const auto& cord : chiral::v_order)
      {
        // Iterate over isospin channels.
        for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
          {
            // Iterate over sectors.
            for (std::size_t sector_index = 0; sector_index < sectors[T0].size(); ++sector_index)
              {
                std::size_t bra_subspace_index = sectors[T0].GetSector(sector_index).bra_subspace_index();
                const basis::RelativeSubspaceLSJT& bra_subspace = sectors[T0].GetSector(sector_index).bra_subspace();
                std::size_t ket_subspace_index = sectors[T0].GetSector(sector_index).ket_subspace_index();
                const basis::RelativeSubspaceLSJT& ket_subspace = sectors[T0].GetSector(sector_index).ket_subspace();

                for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
                  {
                    const basis::RelativeStateLSJT bra_state(bra_subspace, bra_index);
                    for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ ket_index)
                      {
                        const basis::RelativeStateLSJT ket_state(ket_subspace, ket_index);
                        util::OscillatorParameter b(hw);
                        auto rme = op->ReducedMatrixElement(cord, bra_state, ket_state, b,
                                                            regularize, regulator, T0, Abody);
                        temp_matrices[T0][sector_index](bra_index, ket_index) = rme;
                        matrices[T0][sector_index](bra_index, ket_index) += rme;
                      }
                  }

              }
          }

        // Write the contribution at each order.
        std::string cord_str = chiral::reverse_m_order[cord];
        std::string order_file = fmt::format("{}_{:d}n_rel_{}_Nmax{:d}_Jmax{:d}_hw{:.1f}",
                                             name, Abody, cord_str, Nmax, Jmax, hw);
        if(regularize)
          order_file += fmt::format("_regulator_{:.1f}", regulator);
        order_file += fmt::format("_{}.dat", time_str);
        basis::WriteRelativeOperatorLSJT(order_file, space, labels, sectors, temp_matrices, true);

        if (cord == ord)
          break;
      }
    // Write the cumulative matrix element.
    std::string cumulative_file = fmt::format("{}_{:d}n_rel_{}_cumulative_Nmax_{:d}_Jmax_{:d}_hw_{:.1f}",
                                              name, Abody, order, Nmax, Jmax, hw);
    if(regularize)
      cumulative_file += fmt::format("_regulator_{:.1f}", regulator);
    cumulative_file += fmt::format("_{}.dat", time_str);
    basis::WriteRelativeOperatorLSJT(cumulative_file, space, labels, sectors, matrices, true);
  }


  void WriteRelativeCMFiles(const std::string& name,
                            const std::string& order,
                            const std::size_t& Abody,
                            const double& hw,
                            const std::size_t& Nmax,
                            const std::size_t& T0_min,
                            const std::size_t& T0_max,
                            const bool& regularize,
                            const double& regulator)
  {
    // Create operator.
    fmt::print("Generating {} matrix elements...\n", name);
    auto op = chiral::Operator::make(name);
    chiral::Order ord;
    // try
    //   {
    //     ord = chiral::m_order.at(order);
    //   }
    // catch(std::out_of_range& oor)
    //   {
    //     std::cerr << "Chiral order must be in [lo, n4lo].\n";
    //     std::terminate();
    //   }
    ord = chiral::m_order.at(order);

    // Set up operator and file header parameters.
    fmt::print("Beginning RelativeCMLSJT operator basis setup...\n");
    basis::OperatorLabelsJT labels(op->J0(), op->G0(), T0_min, T0_max,
                                   basis::SymmetryPhaseMode::kHermitian);
    basis::RelativeCMOperatorParametersLSJT params(labels, Nmax);

    // Set up relative space.
    basis::RelativeCMSpaceLSJT space(Nmax);

    // Set up operator containers.
    // These arrays store the information for each isospin channel.
    std::array<basis::RelativeCMSectorsLSJT, 3> sectors;
    std::array<basis::OperatorBlocks<double>, 3> matrices;

    // Populate containers.
    for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
      {
        sectors[T0] = basis::RelativeCMSectorsLSJT(space, op->J0(), T0, op->G0());
        basis::SetOperatorToZero(sectors[T0], matrices[T0]);
      }


    // Print truncation information.
    fmt::print("Truncation: Nmax {:d} T0_max {:d}\n", Nmax, T0_max);
    for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
      {
        fmt::print(" T0 {:d} Upper triangular entries {:d} Allocated entries {:d}\n",
                   T0, basis::UpperTriangularEntries(sectors[T0]),
                   basis::AllocatedEntries(matrices[T0]));
      }

    // Time information.
    std::string time_str = GetTimeInfo();

    // Temporary matrices for storing contributions at each order.
    std::array<basis::OperatorBlocks<double>, 3> temp_matrices = matrices;

    // Iterate over chiral orders.
    for (const auto& cord : chiral::v_order)
      {
        // Iterate over isospin channels.
        for (std::size_t T0 = T0_min; T0 <= T0_max; ++T0)
          {
            // Iterate over sectors.
            for (std::size_t sector_index = 0; sector_index < sectors[T0].size(); ++sector_index)
              {
                std::size_t bra_subspace_index = sectors[T0].GetSector(sector_index).bra_subspace_index();
                const basis::RelativeCMSubspaceLSJT& bra_subspace = sectors[T0].GetSector(sector_index).bra_subspace();
                std::size_t ket_subspace_index = sectors[T0].GetSector(sector_index).ket_subspace_index();
                const basis::RelativeCMSubspaceLSJT& ket_subspace = sectors[T0].GetSector(sector_index).ket_subspace();

                for (std::size_t bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
                  {
                    const basis::RelativeCMStateLSJT bra_state(bra_subspace, bra_index);
                    for (std::size_t ket_index = 0; ket_index < ket_subspace.size(); ++ ket_index)
                      {
                        const basis::RelativeCMStateLSJT ket_state(ket_subspace, ket_index);
                        util::OscillatorParameter b(hw);
                        auto rme = op->ReducedMatrixElement(cord, bra_state, ket_state, b,
                                                            regularize, regulator, T0, Abody);
                        temp_matrices[T0][sector_index](bra_index, ket_index) = rme;
                        matrices[T0][sector_index](bra_index, ket_index) += rme;
                      }
                  }

              }

          }
        // Write the contribution at each order.
        std::string cord_str = chiral::reverse_m_order[cord];
        std::string order_file = fmt::format("{}_{:d}n_relcm_{}_Nmax{:d}_hw{:.1f}",
                                             name, Abody, cord_str, Nmax, hw);
        if(regularize)
          order_file += fmt::format("_regulator_{:.1f}", regulator);
        order_file += fmt::format("_{}.dat", time_str);
        basis::WriteRelativeCMOperatorLSJT(order_file, space, labels, sectors, temp_matrices, true);

        if (cord == ord)
          break;
      }
    // Write the cumulative matrix element.
    std::string cumulative_file = fmt::format("{}_{:d}n_relcm_{}_cumulative_Nmax_{:d}_hw_{:.1f}",
                                              name, Abody, order, Nmax, hw);
    if(regularize)
      cumulative_file += fmt::format("_regulator_{:.1f}", regulator);
    cumulative_file += fmt::format("_{}.dat", time_str);
    basis::WriteRelativeCMOperatorLSJT(cumulative_file, space, labels, sectors, matrices, true);
  }
}

#endif
