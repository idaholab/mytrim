#include "../include/TrimVacCount.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimVacCount::TrimVacCount(SimconfType * simconf, SampleBase * sample) :
    TrimBase(simconf, sample),
    _vac_bin()
{
}

void
TrimVacCount::vacancyCreation()
{
  _simconf->vacancies_created++;

  int x = _recoil->_pos(0);
  if (x < 0) return;
  if (x >= int(_vac_bin.size()))
    _vac_bin.resize(x + 1, 0);
  _vac_bin[x]++;
}

void
TrimVacCount::startOutput()
{
  TrimBase::startOutput();
}

void
TrimVacCount::stopOutput()
{
  if (outputting())
  {
    TrimBase::stopOutput();

    // write out vaccancy histogram
    std::ofstream out((_base_name + "_vac.dat").c_str());
    for (unsigned int x = 0; x < _vac_bin.size(); ++x)
      out << x << ' ' << _vac_bin[x] << '\n';
  }
}
