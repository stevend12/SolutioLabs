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
// Gamma 1D Test Program                                                      //
// (gamma_1d.cpp)                                                             //
//                                                                            //
// Steven Dolly                                                               //
// July 28, 2020                                                              //
//                                                                            //
// This program is a test of the 1D Gamma Index calculation function for the  //
// SolutioCpp library. The program replicates the experiment from one of the  //
// initial gamma index journal articles:                                      //
//                                                                            //
// Low DA, Harms WB, Mutic S, Purdy JA. A technique for the quantitative      //
// evaluation of dose distributions. Med Phys. 1998;25(5):656-661.            //
// doi:10.1118/1.598248                                                       //
//                                                                            //
// Result data can be viewed using Octave/MATLAB and the gamma_1d.m script,   //
// located on the Plots directory.                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// C++ headers
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Solutio library headers
#include <SolutioCpp/Therapy/GammaIndex.hpp>
#include <SolutioCpp/Utilities/DataInterpolation.hpp>

int main()
{
  std::cout << "Gamma Index Calculation (1D): Comparison\n";

  ///////////
  // Input //
  ///////////

  const double field_size = 10.0; // Radiation square field size (cm)
  const double profile_width = 20.0; // Profile width (cm)
  const int num_samples = 256; // Number of points in dose profile (size of vectors)

  const double eta = 1.025; // Dose profile scaling (for test profile)
  const double shift = 0.25; // Dose profile shift (for test profile)

  // Set gamma calculation settings
  solutio::GammaIndexSettings s;
  s.GlobalMax = true;    // Use global max
  s.DoseCriteria = 0.03; // 3 %
  s.DistCriteria = 0.3;  // 3 mm
  s.ResampleRate = 0.01; // Resample 100x
  s.Threshold = 0.1;     // 10 %

  /////////////////////////////////////////////////////
  // Step 1: Create test and reference dose profiles //
  /////////////////////////////////////////////////////

  // Dose profile fitting parameters (from reference article)
  const double a = 0.173;
  const double b1 = 0.456;
  const double b2 = 2.892;
  const double t = 0.01;
  const double x0 = 0.5*field_size;
  double x_shift = x0 + shift;

  // Make x-axis
  std::vector<double> x_axis;
  for(int n = 0; n < num_samples; n++)
  {
    x_axis.push_back(-profile_width/2.0 + profile_width*(double(n)/double(num_samples)));
  }

  // Create dose profiles
  std::vector< std::pair<double,double> > test_profile, ref_profile;
  for(int n = 0; n < num_samples; n++)
  {
    std::pair<double,double> temp_dose;
    temp_dose.first = x_axis[n];
    // Reference dose
    temp_dose.second = t + (1.0-t)*(a*((erf(b1*(x0 - fabs(x_axis[n])))+1.0)/2.0) + (1.0-a)*((erf(b2*(x0 - fabs(x_axis[n])))+1.0)/2.0));
    ref_profile.push_back(temp_dose);
    // Test dose (reference dose with some scaling and shift applied)
    temp_dose.second = eta * (t + (1.0-t)*(a*((erf(b1*(x_shift - fabs(x_axis[n])))+1.0)/2.0) + (1.0-a)*((erf(b2*(x_shift - fabs(x_axis[n])))+1.0)/2.0)));
    test_profile.push_back(temp_dose);
  }

  //////////////////////////////////////////////////
  // Step 2. Calculate dose comparison statistics //
  //////////////////////////////////////////////////

  // 2A. Dose difference
  std::vector<double> dose_difference;
  for(int n = 0; n < num_samples; n++)
  {
    dose_difference.push_back(fabs(test_profile[n].second - ref_profile[n].second));
  }

  // 2B. Distance-to-agreement (DTA)
  // Make resampled dose
  double xp;
  solutio::DoublePairVec resampled_ref;
  std::pair<double,double> temp_resample;
  temp_resample.first = xp = ref_profile[0].first;
  temp_resample.second = ref_profile[0].second;
  resampled_ref.push_back(temp_resample);
  do
  {
    xp += 0.001;
    temp_resample.first = xp;
    temp_resample.second = solutio::LinearInterpolation<double>(ref_profile, xp);
    resampled_ref.push_back(temp_resample);
  } while (xp < ref_profile.back().first-0.001);
  // Calculate DTA
  std::vector<double> dta;
  double dd, min_dist;
  bool found;
  for(int m = 0; m < num_samples; m++)
  {
    found = false;
    for(auto it = resampled_ref.begin(); it != resampled_ref.end(); ++it)
    {
      dd = fabs(test_profile[m].second - it->second);
      if(dd <= 0.001)
      {
        if(!found)
        {
          found = true;
          min_dist = fabs(x_axis[m] - it->first);
        }
        else min_dist = std::min(min_dist, fabs(x_axis[m] - it->first));
      }
    }
    if(found) dta.push_back(min_dist);
    else dta.push_back(10.0);
  }

  // 2C. Gamma index and pass rate
  double gamma_pass_rate;
  std::vector<double> results =
    solutio::CalcGammaIndex(test_profile, ref_profile, s, gamma_pass_rate);

  std::cout << std::setprecision(4);
  std::cout << "Gamma Pass Rate (Initial): " << 100.0*gamma_pass_rate << "%\n";

  ////////////////////////////////////////////////
  // Step 3. Analyze gamma calculation settings //
  ////////////////////////////////////////////////

  // Show without global max
  s.GlobalMax = false;
  std::vector<double> results_local_max =
    solutio::CalcGammaIndex(test_profile, ref_profile, s, gamma_pass_rate);
  std::cout << "Gamma Pass Rate (Local Max.): " << 100.0*gamma_pass_rate << "%\n";

  // Show with stricter dose/distance criteria (2 % / 2 mm)
  s.GlobalMax = true;
  s.DoseCriteria = 0.02;
  s.DistCriteria = 0.2;
  std::vector<double> results_1_1 =
    solutio::CalcGammaIndex(test_profile, ref_profile, s, gamma_pass_rate);
  std::cout << "Gamma Pass Rate (2 %, 2 mm): " << 100.0*gamma_pass_rate << "%\n";

  // Show with various resample rates
  s.DoseCriteria = 0.03;
  s.DistCriteria = 0.3;
  s.ResampleRate = 1.0;
  std::vector<double> results_1x =
    solutio::CalcGammaIndex(test_profile, ref_profile, s, gamma_pass_rate);
  std::cout << "Gamma Pass Rate (Resample 1x): " << 100.0*gamma_pass_rate << "%\n";


  /////////////////////////////////////////////////////////////////
  // Step 4. Print profiles to be plotted in Octave (gamma_1d.m) //
  /////////////////////////////////////////////////////////////////
  std::ofstream fout("gamma_1d.txt");
  for(int n = 0; n < num_samples; n++)
  {
    fout << x_axis[n] << ' ' << test_profile[n].second << ' ' <<
      ref_profile[n].second << ' ' << dose_difference[n] << ' ' <<
      dta[n] << ' ' << results[n] << ' ' << results_local_max[n] << ' ' <<
      results_1_1[n] << ' ' << results_1x[n] << '\n';
  }
  fout.close();

  return 0;
}
