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

#include "../include/TrimVacCount.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimVacCount::TrimVacCount(SimconfType * simconf, SampleBase * sample)
  : ThreadedTrimBase(simconf, sample), _vac_bin()
{
}

void
TrimVacCount::vacancyCreation()
{
  _simconf->vacancies_created++;

  int x = _recoil->_pos(0);
  if (x < 0)
    return;
  if (x >= int(_vac_bin.size()))
    _vac_bin.resize(x + 1, 0);
  _vac_bin[x]++;
}

void
TrimVacCount::replacementCollision()
{
  int x = _recoil->_pos(0);
  if (x < 0)
    return;
  if (x >= int(_repl_bin.size()))
    _repl_bin.resize(x + 1, 0);
  _repl_bin[x]++;
}

void
TrimVacCount::threadJoin(const ThreadedTrimBase & ttb)
{
  const TrimVacCount & tvc = static_cast<const TrimVacCount &>(ttb);

  // sum histogram data
  _vac_bin.resize(std::max(_vac_bin.size(), tvc._vac_bin.size()));
  for (unsigned int x = 0; x < tvc._vac_bin.size(); ++x)
    _vac_bin[x] += tvc._vac_bin[x];
  _repl_bin.resize(std::max(_repl_bin.size(), tvc._repl_bin.size()));
  for (unsigned int x = 0; x < tvc._repl_bin.size(); ++x)
    _repl_bin[x] += tvc._repl_bin[x];
}

void
TrimVacCount::writeOutput()
{
  // write out vaccancy histogram
  const unsigned int size = std::max(_vac_bin.size(), _repl_bin.size());
  _vac_bin.resize(size);
  _repl_bin.resize(size);
  std::ofstream out((_base_name + "_vac.dat").c_str());
  for (unsigned int x = 0; x < size; ++x)
    out << x << ' ' << _vac_bin[x] << ' ' << _repl_bin[x] << '\n';
}
