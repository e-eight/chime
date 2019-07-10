#include <iomanip>
// #include <chrono>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <string>
#include "fmt/format.h"
#include "fmt/time.h"
#include "fmt/chrono.h"
#include "fmt/ostream.h"
#include "basis/lsjt_operator.h"
#include "chiral.h"
#include "constants.h"
#include "cli.h"
#include "operators.h"

int main(int argc, char** argv)
{
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////// String to Order map ////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  std::vector<std::string> order_name;
  std::unordered_map<std::string, chiral::Order> order_map;

  order_map["lo"] = chiral::Order::lo;
  order_name.push_back("lo");
  order_map["nlo"] = chiral::Order::nlo;
  order_name.push_back("nlo");
  order_map["n2lo"] = chiral::Order::n2lo;
  order_name.push_back("n2lo");
  order_map["n3lo"] = chiral::Order::n3lo;
  order_name.push_back("n3lo");
  order_map["n4lo"] = chiral::Order::n4lo;
  order_name.push_back("n4lo");

  //////////////////////////////////////////////////////////////////////////////
  /////////////// Application description and input parameters /////////////////
  //////////////////////////////////////////////////////////////////////////////

  CLI::App app("Generates CEFT reduced matrix elements in HO basis.");

  // flags
  std::string name = "rsq";
  std::string order = "lo";
  double hw = 10;
  int Nmax = 200;
  int Jmax = 1;
  int T0_min = 0;
  int T0_max = 0;
  double regulator = 0.9; // (LENPIC regulator in fm)

  app.add_option("-n,--name", name, "Name of operator.");
  app.add_option("-o,--order", order, "Chiral order of operator.");
  app.add_option("-E,--hw", hw, "Oscillator energy of basis.");
  app.add_option("-N,--Nmax", Nmax, "Nmax truncation of basis.");
  app.add_option("-J,--Jmax", Jmax, "Jmax truncation of basis.");
  app.add_option("-t,--T0_min", T0_min, "Minimum isospin component of operator.");
  app.add_option("-T,--T0_max", T0_max, "Maximum isospin component of operator.");
  app.add_option("-R,--regulator", regulator, "Value of LENPIC regulator, in fermi.");

  app.set_config("-c,--config");

  // Parse input
  try
    {
      app.parse(argc, argv);
    }
  catch (const CLI::ParseError &e)
    {
      return app.exit(e);
    }

  //////////////////////////////////////////////////////////////////////////////
  /////////////////// Create chiral operator from input ////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  auto op = chiral::Operator::make(name);

  // Print Header
  fmt::print("Generating {} matrix elements...\n", name);

  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////// Create relative basis //////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  fmt::print("Beginning RelativeLSJT operator basis setup...\n");

  // Set operator and file header parameters
  basis::OperatorLabelsJT op_labels(op->J0(), op->G0(), T0_min, T0_max,
                                    basis::SymmetryPhaseMode::kHermitian);
  basis::RelativeOperatorParametersLSJT op_params(op_labels, Nmax, Jmax);

  // Set up relative space
  basis::RelativeSpaceLSJT space(op_params.Nmax, op_params.Jmax);

  // Set up operator containers
  // These are arrays to store T0 = 0/1/2 components.
  std::array<basis::RelativeSectorsLSJT, 3> sectors;
  std::array<basis::OperatorBlocks<double>, 3> matrices;

  // Populate operator containers
  basis::ConstructZeroOperatorRelativeLSJT(op_params, space, sectors, matrices);


  //////////////////////////////////////////////////////////////////////////////
  /////////////////////// Generate matrix elements /////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  // Operator diagnostics
  fmt::print("Truncation: Nmax {0} Jmax {1} T0_max {2}\n",
             op_params.Nmax, op_params.Jmax, op_params.T0_max);

  fmt::print("Matrix elements:");
  for (auto T0 = op_params.T0_min; T0 <= op_params.T0_max; ++T0)
    fmt::print(" {}\n", basis::UpperTriangularEntries(sectors[T0]));
  fmt::print("Allocated:");
  for (auto T0 = op_params.T0_min; T0 <= op_params.T0_max; ++T0)
    fmt::print(" {}\n", basis::AllocatedEntries(matrices[T0]));

  // Populate matrix elements

  // Oscillator length scale
  const auto hbarc2 = constants::hbarc * constants::hbarc;
  const auto osc_b = std::sqrt(hbarc2 / constants::reduced_nucleon_mass_MeV / hw);
  // Oscillator energy string for filename
  std::string hw_str = fmt::format("{:.1f}", hw);

  // Get time for output filename
  auto now = std::chrono::system_clock::now();
  auto time = std::chrono::system_clock::to_time_t(now);
  auto time_str = std::to_string(time);

  // Temporary matrices for storing matrix elements of each order
  std::array<basis::OperatorBlocks<double>, 3> temp_matrices = matrices;

  // Iterate over chiral orders
  for (auto current_order : order_name)
    {
      // Iterate over isospin
      for (auto T0 = op_params.T0_min; T0 <= op_params.T0_max; ++T0)
        {
          // Iterate over sectors
          for (int sector_index = 0; sector_index < sectors[T0].size(); ++sector_index)
            {
              // Get bra and ket subspaces
              const basis::RelativeSectorsLSJT::SectorType& sector =
                sectors[T0].GetSector(sector_index);
              const basis::RelativeSectorsLSJT::SubspaceType& bra_subspace =
                sector.bra_subspace();
              const basis::RelativeSectorsLSJT::SubspaceType& ket_subspace =
                sector.ket_subspace();

              // Get states
              for (int bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
                {
                  const basis::RelativeStateLSJT bra_state(bra_subspace, bra_index);
                  for (int ket_index = 0; ket_index < ket_subspace.size(); ++ket_index)
                    {
                      const basis::RelativeStateLSJT ket_state(ket_subspace, ket_index);

                      // Calculate matrix element
                      auto order_enum = order_map[current_order];

                      temp_matrices[T0][sector_index](bra_index, ket_index) =
                        op->ReducedMatrixElement(order_enum, bra_state, ket_state, osc_b, regulator);

                      matrices[T0][sector_index](bra_index, ket_index) +=
                        op->ReducedMatrixElement(order_enum, bra_state, ket_state, osc_b, regulator);
                    }
                }
            }
        }
      // Write the contributions at each order
      std::string order_file =
        name + "_2b_rel_" + current_order + "_N" + std::to_string(op_params.Nmax)
          + "_J" + std::to_string(op_params.Jmax) + "_hw" + hw_str
        + "_" + time_str + ".dat";
      basis::WriteRelativeOperatorLSJT(order_file, space, op_labels,
                                       sectors, temp_matrices, true);

      // Break if required order has been reached
      if (current_order == order)
        break;
    }

  // Write cumulative contribution
  std::string cumulative_file =
    name + "_2b_rel_" + order + "_cumulative"
    + "_N" + std::to_string(op_params.Nmax)
    + "_J" + std::to_string(op_params.Jmax)
    + "_hw" + hw_str + "_" + time_str + ".dat";
  basis::WriteRelativeOperatorLSJT(cumulative_file, space, op_labels,
                                   sectors, matrices, true);
}
