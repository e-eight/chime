/*******************************************************************************
 constants.h

 Defines the various physical and mathematical constants needed for the matrix
 element calculations.

 Language: C++11
 Soham Pal
 Iowa State University
*******************************************************************************/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

namespace chime {
namespace constants {

constexpr double pi = 3.141592653589793238462643383279502884;

// Conversion Factor
constexpr double hbarc = 197.3269602;       // (in MeV fm)
constexpr double hbarc_GeV = 0.1973269602;  // (in GeV fm)

// Masses
// The units are either MeV or fm^{-1}.
constexpr double charged_pion_mass_MeV = 139.57061;
constexpr double neutral_pion_mass_MeV = 134.9770;
constexpr double pion_mass_MeV =
    (2 * charged_pion_mass_MeV + neutral_pion_mass_MeV) / 3.;
constexpr double pion_mass_fm = pion_mass_MeV / hbarc;
constexpr double proton_mass_MeV = 938.27208816;
constexpr double proton_mass_fm = proton_mass_MeV / hbarc;
constexpr double neutron_mass_MeV = 939.56542052;
constexpr double neutron_mass_fm = neutron_mass_MeV / hbarc;
constexpr double nucleon_mass_MeV = ((proton_mass_MeV + neutron_mass_MeV) / 2);
constexpr double nucleon_mass_fm = nucleon_mass_MeV / hbarc;
constexpr double nuclear_magneton_MeV =
    1.0 / (nucleon_mass_MeV * 2);                                    // (e = 1)
constexpr double nuclear_magneton_fm = 1.0 / (nucleon_mass_fm * 2);  // (e = 1)

// Low energy constants
constexpr double isoscalar_nucleon_charge_radius_sq_fm = 0.603729;  // (in fm^2)
constexpr double isoscalar_nucleon_magnetic_moment = 0.8798046315;
constexpr double isovector_nucleon_magnetic_moment = 4.7058900815;
constexpr double gA = 1.29;
constexpr double pion_decay_constant_MeV = 92.4;
constexpr double pion_decay_constant_fm =
    pion_decay_constant_MeV / hbarc;  // (in fm^{-1})
constexpr double d9_GeV = -0.06;  // (in GeV^{-2}, from Gasparyan, Lutz, 2010)
constexpr double d9_fm = d9_GeV * hbarc_GeV * hbarc_GeV;  // (in fm^2)

// constexpr double L2_GeV = 0.188; // (in GeV^{-4})
constexpr double L2_fm = -0.149;  // (in fm^4, from Chen, Rupak, Savage, 1999)
constexpr double c3_GeV = -4.69;  // (in GeV^{-1})
constexpr double c3_fm = c3_GeV * hbarc_GeV;
constexpr double c4_GeV = 3.40;  // (in GeV^{-1})
constexpr double c4_fm = c4_GeV * hbarc_GeV;
constexpr double cD = 2.1;  // (Recommended value from Epelbaum et.al. 2018)

}  // namespace constants
}  // namespace chime

#endif
