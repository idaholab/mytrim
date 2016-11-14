#include "../include/TrimVacEnergyCount.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimVacEnergyCount::TrimVacEnergyCount(SimconfType * simconf, SampleBase * sample) :
    ThreadedTrimBase(simconf, sample),
    _evac_bin()
{
}

void
TrimVacEnergyCount::vacancyCreation()
{
  _simconf->vacancies_created++;

  // position bin
  int x = _recoil->_pos(0);
  if (x < 0) return;

  // natural logarithm of the recoil energy
  int E = std::log(_recoil->_E);
  E = E < 0 ? 0 : E;

  // dynamically resize the arrays
  if (E >= int(_evac_bin.size()))
    _evac_bin.resize(E + 1);
  if (x >= int(_evac_bin[E].size()))
    _evac_bin[E].resize(x + 1, 0);

  // increment counter
  _evac_bin[E][x]++;
}

void
TrimVacEnergyCount::threadJoin(const ThreadedTrimBase & ttb)
{
  const TrimVacEnergyCount & tvec = static_cast<const TrimVacEnergyCount &>(ttb);

  // sum histogram data
  _evac_bin.resize(std::max(_evac_bin.size(), tvec._evac_bin.size()));
  for (unsigned int E = 0; E < tvec._evac_bin.size(); ++E)
  {
    _evac_bin[E].resize(std::max(_evac_bin[E].size(), tvec._evac_bin[E].size()));
    for (unsigned int x = 0; x < tvec._evac_bin[E].size(); ++x)
      _evac_bin[E][x] += tvec._evac_bin[E][x];
  }
}

void
TrimVacEnergyCount::writeOutput()
{
  // write out vaccancy histogram
  std::ofstream out((_base_name + "_evac.dat").c_str());
  for (unsigned int E = 0; E < _evac_bin.size(); ++E)
  {
    for (unsigned int x = 0; x < _evac_bin[E].size(); ++x)
      out << E << ' ' << x << ' ' << _evac_bin[E][x] << '\n';
    out << '\n';
  }
}
