/******************************************************************************\
MIT License

Copyright (c) 2020 Steven Dolly

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
\******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
// CT HU Calculation Program                                                  //
// (ct_hu.cpp)                                                                //
//                                                                            //
// Steven Dolly                                                               //
// July 28, 2020                                                              //
//                                                                            //
// This program calculates CT HU for various materials common in medical CT   //
// scans. It combines the Tasmip and NistPad functions of the SolutioCpp      //
// library to calculate HU ranges for these materials.                        //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

// Solutio library headers
#include <SolutioCpp/Physics/NistPad.hpp>
#include <SolutioCpp/Imaging/Tasmip.hpp>

int main(int argc, char *argv[])
{
  std::cout << "CT HU Calculation\n\n";

  std::string folder = std::string(argv[1]);

  // Calculate tungsten x-ray spectum using the TASMIP algorithm
  std::vector<double> kVp80 = solutio::Tasmip(80, 0.0, "Aluminum", folder);
  std::vector<double> kVp100 = solutio::Tasmip(100, 0.0, "Aluminum", folder);
  std::vector<double> kVp120 = solutio::Tasmip(120, 0.0, "Aluminum", folder);
  std::vector<double> kVp140 = solutio::Tasmip(140, 0.0, "Aluminum", folder);

  // Calculate mean energies for each spectrum; print to terminal
  double mean_energy[4] = {0.0, 0.0, 0.0, 0.0};
  for(int n = 0; n < 151; n++)
  {
    mean_energy[0] += kVp80[n]*n;
    mean_energy[1] += kVp100[n]*n;
    mean_energy[2] += kVp120[n]*n;
    mean_energy[3] += kVp140[n]*n;
  }
  std::cout << "Mean Energies\n";
  std::cout << "-------------\n";
  std::cout << "80 kVp: " << mean_energy[0] << "\n";
  std::cout << "100 kVp: " << mean_energy[1] << "\n";
  std::cout << "120 kVp: " << mean_energy[2] << "\n";
  std::cout << "140 kVp: " << mean_energy[3] << "\n\n";

  // Load NistPad sample materials
  solutio::NistPad NistAir(folder, "Air");
  solutio::NistPad NistWater(folder, "Water");

  const int num_materials = 9;
  std::string material_names[num_materials] = {"Air", "Water", "Adipose",
    "Bone", "Brain", "Breast", "EyeLens", "Lung", "Muscle"};
  std::vector<solutio::NistPad> materials;

  std::cout << "Material & Density (g/cm^3)\n";
  std::cout << "---------------------------\n";
  for(int n = 0; n < num_materials; n++)
  {
    materials.push_back(solutio::NistPad(folder, material_names[n]));
    if(n == 7) materials[n].ForceDensity(0.25);
    std::cout << material_names[n] << ": " << materials[n].GetDensity() << ", "
      << materials[n].PowerLawEffectiveZ(3.0) << '\n';
  }
  std::cout << '\n';

  // Calculate CT numbers for each material and spectrum
  std::cout << "CT HU By Material\n";
  std::cout << "-----------------\n";
  int ct_hu[num_materials][2];
  for(int m = 0; m < num_materials; m++)
  {
    for(int n = 0; n < 4; n++)
    {
      int hu = round( 1000.0 *
        ((materials[m].LinearAttenuation(mean_energy[n]) -
        NistWater.LinearAttenuation(mean_energy[n])) /
        (NistWater.LinearAttenuation(mean_energy[n]) -
        NistAir.LinearAttenuation(mean_energy[n]))));
      if(n == 0)
      {
        ct_hu[m][0] = ct_hu[m][1] = hu;
      }
      else
      {
        if(hu < ct_hu[m][0]) ct_hu[m][0] = hu;
        if(hu > ct_hu[m][1]) ct_hu[m][1] = hu;
      }
    }
    std::cout << material_names[m] << ": " << ct_hu[m][0] << " -> " <<
      ct_hu[m][1] << '\n';
  }

  // Print results to file for plotting using Octave/MATLAB
  std::ofstream fout("ct_hu.txt");
  for(int n = 0; n < num_materials; n++)
  {
    fout << material_names[n] << "," << materials[n].GetDensity() << "," <<
      materials[n].PowerLawEffectiveZ(3.0) << "," << ct_hu[n][0] << "," <<
      ct_hu[n][1] << '\n';
  }
  fout.close();

  return 0;
}
