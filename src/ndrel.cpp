#include <iostream>
#include <iomanip>
#include "basis/lsjt_scheme.h"
#include "basis/lsjt_operator.h"
#include "chiral.h"
#include "operator_rc_sq.h"
#include "op_file_io.h"
#include "constants.h"

int main()
{
    // From input to parameters for matrix generation.
    auto params = chiral::input_to_params();
    auto name = params.name, order = params.order;
    auto op = std::move(params.op);
    auto J0_min = params.J0_min, J0_max = params.J0_max;
    auto T0_min = params.T0_min, T0_max = params.T0_max;
    auto Nmax = params.Nmax, Jmax = params.Jmax;
    auto hw = params.hw;
    
    // Oscillator constant and hw string
    const double osc_b = std::sqrt(constants::HBARC * constants::HBARC
                                   / constants::RED_NUCLEON_MASS / hw);
    std::ostringstream hw_str;
    hw_str << hw;

    // Print header
    std::cout << "\n";
    std::cout << "Now generating operator: " << name << "\n";

    // Set up zero operator for relative format
    std::cout << "Beginning RelativeLSJT operator basis setup..." << "\n";
    basis::RelativeSpaceLSJT space(Nmax, Jmax);
    basis::OperatorLabelsJT
        operator_labels(op->J0,
                        op->G0,
                        T0_min,
                        T0_max,
                        basis::SymmetryPhaseMode::kHermitian);
    basis::RelativeOperatorParametersLSJT operator_parameters(operator_labels,
                                                              Nmax,
                                                              Jmax);
    std::array<basis::RelativeSectorsLSJT, 3> sectors;
    std::array<basis::OperatorBlocks<double>, 3> matrices;
    basis::ConstructZeroOperatorRelativeLSJT(operator_parameters,
                                             space,
                                             sectors,
                                             matrices);
    // Operator diagonostics
    std::cout << " Truncation: "
              << " Nmax " << Nmax
              << " Jmax " << Jmax
              << " T0_max " << T0_max
              << "\n";
    std::cout << " Matrix elements: ";
    for (auto T0 = T0_min; T0 <= T0_max; ++T0)
        std::cout <<  basis::UpperTriangularEntries(sectors[T0]);
    std::cout << "\n";
    std::cout << " Allocated: ";
    for (auto T0 = T0_min; T0 <= T0_max; ++T0)
        std::cout << basis::AllocatedEntries(matrices[T0]);
    std::cout << "\n";

    // Populate reduced matrix elements
    // Iterate over T0
    for (int t0_index = 0; t0_index < 2; ++t0_index)
    {
        for (int sector_index = 0;
             sector_index < sectors[t0_index].size(); ++sector_index)
        {
            const typename basis::RelativeSectorsLSJT::SectorType& sector =
                sectors[t0_index].GetSector(sector_index);
            const typename basis::RelativeSectorsLSJT::SubspaceType&
                bra_subspace = sector.bra_subspace();
            const typename basis::RelativeSectorsLSJT::SubspaceType&
                ket_subspace = sector.ket_subspace();

            // Get states
            for (int bra_index = 0; bra_index < bra_subspace.size(); ++bra_index)
            {
                const basis::RelativeStateLSJT bra(bra_subspace, bra_index);
                for (int ket_index = 0;
                     ket_index < ket_subspace.size(); ++ket_index)
                {
                    const basis::RelativeStateLSJT ket(ket_subspace, ket_index);

                    // Calculate the reduced matrix element
                    double reduced_matrix_element =
                        op->calculate_rme(bra, ket, osc_b);

                    // Save reduced matrix element
                    matrices[t0_index][sector_index](bra_index, ket_index)
                        += reduced_matrix_element;
                }
            }
        }
    }

    // Construct output filename
    std::string filename =
        name + "_tbrel_" + order + "_N" + std::to_string(Nmax)
        + "_J" + std::to_string(Jmax) + "_hw" + hw_str.str() + ".dat";

    // Write to output file
    basis::WriteRelativeOperatorLSJT(filename,
                                     space,
                                     operator_labels,
                                     sectors,
                                     matrices,
                                     true);
}
