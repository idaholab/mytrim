/*
MyTRIM - a three dimensional binary collision Monte Carlo library.
Copyright (C) 2008-2018  Daniel Schwen <daniel@schwen.de>

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA
*/

#include "../include/TrimRange.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimRange::TrimRange(SimconfType * simconf, SampleBase * sample)
  : ThreadedTrimBase(simconf, sample), _range(112) // reserve space for all physical Z
{
}

void
TrimRange::vacancyCreation()
{
  const Real Ed = _element->_Edisp;
  const Real ed = 0.0115 * std::pow(_recoil->_Z, -7.0 / 3.0) * _recoil->_E;
  const Real kd = 0.1337 * std::pow(_recoil->_Z, 2.0 / 3.0) / std::sqrt(_recoil->_m);
  const Real g = 3.4008 * std::pow(ed, 1.0 / 6.0) + 0.40244 * std::pow(ed, 3.0 / 4.0) + ed;
  const Real Ev = _recoil->_E / (1.0 + kd * g);

  // Quick calculation of Damage
  if (Ev < Ed)
    return;
  if (Ev >= Ed / 0.4)
    _simconf->vacancies_created += Ev * 0.4 / Ed;
  else
    _simconf->vacancies_created++;
}

void
TrimRange::dissipateRecoilEnergy()
{
  // store the x-component of the stopped ion position
  _range[_recoil->_Z].push_back(_recoil->_pos(0));
}

void
TrimRange::threadJoin(const ThreadedTrimBase & ttb)
{
  const TrimRange & tr = static_cast<const TrimRange &>(ttb);

  // combine range data
  for (unsigned int Z = 0; Z < _range.size(); ++Z)
    _range[Z].insert(_range[Z].end(), tr._range[Z].begin(), tr._range[Z].end());
}

void
TrimRange::writeOutput()
{
  // determine histogram range and bins
  unsigned max_samples = 0;
  Real x_min = 0.0;
  Real x_max = 0.0;
  for (unsigned int Z = 0; Z < _range.size(); ++Z)
  {
    if (_range[Z].size() > max_samples)
      max_samples = _range[Z].size();

    for (auto x : _range[Z])
    {
      x_min = std::min(x_min, x);
      x_max = std::max(x_max, x);
    }
  }

  // heuristic for bin size
  const Real xwidth = x_max - x_min;
  const Real bin = std::min(xwidth * 100.0 / max_samples, xwidth / 10.0);
  const unsigned int nbin = xwidth / bin + 1;

  // build all histograms
  std::vector<std::pair<int, std::vector<unsigned int>>> histograms;
  for (unsigned int Z = 0; Z < _range.size(); ++Z)
  {
    if (_range[Z].empty())
      continue;

    std::vector<unsigned int> histogram(nbin);
    for (auto x : _range[Z])
      histogram[std::floor((x - x_min) / bin)]++;

    histograms.push_back(std::make_pair(Z, histogram));
  }

  // open range histogram for output
  std::ofstream out((_base_name + "_ranges.dat").c_str());

  // write column header
  out << "#x";
  for (auto & histogram : histograms)
    out << " Z" << histogram.first;
  out << '\n';

  // write data
  for (unsigned int n = 0; n < nbin; ++n)
  {
    out << (n * bin + x_min);
    for (auto & histogram : histograms)
      out << ' ' << histogram.second[n];
    out << '\n';
  }
}
