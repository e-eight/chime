#include <iomanip>
#include "io_new.h"
#include "CLI/CLI.hpp"

int main(int argc, char** argv)
{
  //////////////////////////////////////////////////////////////////////////////
  /////////////// Application description and input parameters /////////////////
  //////////////////////////////////////////////////////////////////////////////

  CLI::App app("Generates CEFT reduced matrix elements in HO basis.");

  // flags
  std::string name = "m1";
  std::string order = "lo";
  std::size_t Abody = 1;
  bool has_cm = false;
  double hw = 10;
  int Nmax = 200;
  int Jmax = 1;
  int T0_min = 0;
  int T0_max = 0;
  bool regularize = false;
  double regulator = 0.9; // (LENPIC regulator in fm)


  app.add_option("-n,--name", name, "Name of operator.");
  app.add_option("-o,--order", order, "Chiral order of operator.");
  app.add_option("-A,--Abody", Abody, "One body or two body chiral operator.");
  app.add_flag("-C,--has_cm", has_cm, "Does the operator depend on the CM?");
  app.add_option("-E,--hw", hw, "Oscillator energy of basis.");
  app.add_option("-N,--Nmax", Nmax, "Nmax truncation of basis.");
  app.add_option("-J,--Jmax", Jmax, "Jmax truncation of basis.");
  app.add_option("-t,--T0_min", T0_min, "Minimum isospin component of operator.");
  app.add_option("-T,--T0_max", T0_max, "Maximum isospin component of operator.");
  app.add_flag("-r,--regularize", regularize, "Do you want the regulated results?");
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
  //////////////////////////// Write RMEs to file //////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  if (has_cm)
    io::WriteRelativeCMFiles(name, order, Abody, hw, Nmax, T0_min, T0_max, regularize, regulator);
  else
    io::WriteRelativeFiles(name, order, Abody, hw, Nmax, Jmax, T0_min, T0_max, regularize, regulator);
}
