#include "../include/TrimVacEnergyCount.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimVacEnergyCount::TrimVacEnergyCount(SimconfType * simconf, SampleBase * sample) :
    TrimBase(simconf, sample),
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
TrimVacEnergyCount::startOutput()
{
  TrimBase::startOutput();
}

void
TrimVacEnergyCount::stopOutput()
{
  if (outputting())
  {
    TrimBase::stopOutput();

    // write out vaccancy histogram
    std::ofstream out((_base_name + "_evac.dat").c_str());
    for (unsigned int E = 0; E < _evac_bin.size(); ++E)
    {
      for (unsigned int x = 0; x < _evac_bin[E].size(); ++x)
        out << E << ' ' << x << ' ' << _evac_bin[E][x] << '\n';
      out << '\n';
    }
  }
}
