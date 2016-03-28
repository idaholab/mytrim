#include "../include/TrimVacEnergyCount.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimVacEnergyCount::TrimVacEnergyCount(SimconfType * simconf, SampleBase * sample) :
    TrimBase(simconf, sample),
    _vac_bin()
{
}

void
TrimVacEnergyCount::vacancyCreation()
{
  _simconf->vacancies_created++;

  // position bin
  int x = _recoil->_pos(0);
  if (x < 0) return;

  // decimal logarithm of the recoil energy
  int E = std::log10(_recoil->_E);
  E = E < 0 ? 0 : E;

  // dynamically resize the arrays
  if (E >= _vac_bin.size())
    _vac_bin.resize(E + 1);
  if (x >= _vac_bin[E].size())
    _vac_bin[E].resize(x + 1, 0);

  // increment counter
  _vac_bin[E][x]++;
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
    for (unsigned int E = 0; E < _vac_bin.size(); ++E)
    {
      for (unsigned int x = 0; x < _vac_bin[E].size(); ++x)
        out << E << ' ' << x << ' ' << _vac_bin[E][x] << '\n';
      out << '\n';
    }
  }
}
