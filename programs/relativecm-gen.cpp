/*******************************************************************************
 relativecm-gen.cpp

 Generates relativecm matrix elements for defined operators. Currently only the
 magnetic moment operator has been defined.

 Standard input:
   J0 g0 T0_min T0_max
   Nmax hw
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
#include "mcutils/parsing.h"
#include "relativecm_rme.h"

// Relative-cm input parameters.
struct InputParameters {
  InputParameters() {}
  InputParameters(std::string input_filename);

  basis::RelativeCMOperatorParametersLSJT op_params;
  double hbomega;
  double R;
  std::string op_name;
  std::string order;
  int abody;
  std::string target_filename;
};

// Read parameters.
InputParameters::InputParameters(std::string input_filename) {
  std::string line;
  int line_count = 0;

  std::ifstream input_file(input_filename);

  if (input_file.is_open()) {
    std::cout << "Reading input parameters... " << "\n";
    while (input_file.good()) {
      // Line 1. Get tensor properties.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> op_params.J0
                    >> op_params.g0
                    >> op_params.T0_min
                    >> op_params.T0_max;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }

      // Line 2. Get operator basis parameters.
      {
        ++line_count;
        std::getline(input_file, line);
        std::istringstream line_stream(line);
        line_stream >> op_params.Nmax >> hbomega;
        mcutils::ParsingCheck(line_stream, line_count, line);
      }

      // Set miscellaneous op_params fields.
      op_params.symmetry_phase_mode = basis::SymmetryPhaseMode::kHermitian;

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
        line_stream >> op_name >> order >> abody;
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
  } else {
    std::runtime_error("File not available.");
  }
}


// Populate operator.
void PopulateOperator(const InputParameters& input_params,
                      basis::RelativeCMSpaceLSJT& relcm_space,
                      std::array<basis::RelativeCMSectorsLSJT, 3>& relcm_sectors,
                      std::array<basis::OperatorBlocks<double>, 3>& relcm_matrices) {

  std::cout << "Populating operator..." << "\n";

  // Set up operator parameters.
  const basis::RelativeCMOperatorParametersLSJT op_params = input_params.op_params;

  // Set up relative-cm space.
  relcm_space = basis::RelativeCMSpaceLSJT(op_params.Nmax);

  // Populate operator containers.
  if (input_params.op_name == "mm") {
    if (input_params.order == "nlo") {
      if (input_params.abody == 2) {
        relcm::ConstructMu2nNLOOperator(op_params,
                                        relcm_space,
                                        relcm_sectors,
                                        relcm_matrices,
                                        input_params.hbomega,
                                        input_params.R);
      }
    }
  }
}

int main()
{
  // Read parameters.
  InputParameters input_params("relcm.in");
  std::cout << "  Operator "
            << input_params.op_name << " "
            << input_params.order << " "
            << input_params.abody << "n\n";
  std::cout << "  Nmax "
            << input_params.op_params.Nmax
            << " hw "
            << input_params.hbomega
            << " R " << input_params.R << "\n";
  std::cout << "  T0_min "
            << input_params.op_params.T0_min
            << " T0_max "
            << input_params.op_params.T0_max << "\n";

  // Set up operator.
  basis::RelativeCMSpaceLSJT relcm_space;
  std::array<basis::RelativeCMSectorsLSJT, 3> relcm_sectors;
  std::array<basis::OperatorBlocks<double>, 3> relcm_matrices;
  PopulateOperator(input_params, relcm_space, relcm_sectors, relcm_matrices);

  // Write operator.
  basis::WriteRelativeCMOperatorLSJT(input_params.target_filename,
                                     relcm_space,
                                     input_params.op_params,
                                     relcm_sectors,
                                     relcm_matrices,
                                     true);
}
