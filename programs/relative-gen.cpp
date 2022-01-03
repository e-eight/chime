/*******************************************************************************
 relative-gen.cpp

 Generates relative matrix elements for defined operators.

 Input (relative.in):
   J0 g0 T0_min T0_max
   Nmax Jmax hw
   R
   operator_name chiral_order one_or_two_body
   output_filename

 The operator_name may be:

   mm
     Magnetic moment operator

 The chiral_order may be:

   lo, nlo, n2lo

 The one_or_two_body may be:

   1
     One body current

   2
     Two body current

 Language: C++11
 Soham Pal
 Iowa State University
*******************************************************************************/

#include <fstream>

#include "chime.h"
#include "mcutils/parsing.h"
#include "relative_rme.h"

// Input parameters for relative operators.
struct InputParameters {
  InputParameters();
  InputParameters(std::string input_filename);

  basis::RelativeOperatorParametersLSJT basis_params;
  double hbomega;
  double R;
  std::string op_name;
  std::string op_order;
  std::size_t op_abody;
  std::string target_filename;
};

InputParameters::InputParameters(std::string input_filename)
{
  std::string line;
  int line_count = 0;

  std::ifstream input_file(input_filename);

  if (input_file.is_open()) {
    std::cout << "Reading input parameters...\n";
    while (input_file.good()) {
      // Line 1. Get tensor properties.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> basis_params.J0 >> basis_params.g0 >> basis_params.T0_min
            >> basis_params.T0_max;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }

      // Line 2. Get operator basis parameters.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> basis_params.Nmax >> basis_params.Jmax >> hbomega;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }

      // Set miscellaneous labels fields.
      basis_params.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

      // Line 3. Get LENPIC semilocal coordinate space regulator parameter.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> R;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }

      // Line 4. Get operator choice.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> op_name >> op_order >> op_abody;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }

      // Line 5. Get output filename.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> target_filename;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }
    }
  }
  else {
    std::runtime_error("File not available.");
  }
}

// Populate operator.
void PopulateOperator(
    const InputParameters &input_params, basis::RelativeSpaceLSJT &rel_space,
    std::array<basis::RelativeSectorsLSJT, 3> &rel_sectors,
    std::array<basis::OperatorBlocks<double>, 3> &rel_matrices)
{
  std::cout << "Populating operator...\n";

  // Set up relative-cm space.
  rel_space = basis::RelativeSpaceLSJT(input_params.basis_params.Nmax,
                                       input_params.basis_params.Jmax);

  // Populate operator containers.
  if (input_params.op_name == "mm") {
    if (input_params.op_order == "nlo") {
      if (input_params.op_abody == 2) {
        chime::relative::ConstructMu2nNLOOperator(
            input_params.basis_params, rel_space, rel_sectors, rel_matrices,
            input_params.hbomega, input_params.R);
      }
    }
  }
}

int main()
{
  // Read parameters.
  InputParameters input_params("relative.in");
  std::cout << "  Operator " << input_params.op_name << " "
            << input_params.op_order << " " << input_params.op_abody << "n\n";
  std::cout << "  Nmax " << input_params.basis_params.Nmax << " Jmax "
            << input_params.basis_params.Jmax << " hw " << input_params.hbomega
            << " R " << input_params.R << "\n";
  std::cout << "  T0_min " << input_params.basis_params.T0_min << " T0_max "
            << input_params.basis_params.T0_max << "\n";

  // Set up operator.
  basis::RelativeSpaceLSJT rel_space;
  std::array<basis::RelativeSectorsLSJT, 3> rel_sectors;
  std::array<basis::OperatorBlocks<double>, 3> rel_matrices;
  PopulateOperator(input_params, rel_space, rel_sectors, rel_matrices);

  // Write operator.
  basis::WriteRelativeOperatorLSJT(input_params.target_filename, rel_space,
                                   input_params.basis_params, rel_sectors,
                                   rel_matrices, true);
}
