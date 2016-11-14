#include "../include/TrimVacCount.h"
#include <fstream>

using namespace MyTRIM_NS;

TrimVacCount::TrimVacCount(SimconfType * simconf, SampleBase * sample) :
    ThreadedTrimBase(simconf, sample),
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
TrimVacCount::threadJoin(const ThreadedTrimBase & ttb)
{
  const TrimVacCount & tvc = static_cast<const TrimVacCount &>(ttb);

  // sum histogram data
  _vac_bin.resize(std::max(_vac_bin.size(), tvc._vac_bin.size()));
  for (unsigned int x = 0; x < tvc._vac_bin.size(); ++x)
    _vac_bin[x] += tvc._vac_bin[x];
}

void
TrimVacCount::writeOutput()
{
  // write out vaccancy histogram
  std::ofstream out((_base_name + "_vac.dat").c_str());
  for (unsigned int x = 0; x < _vac_bin.size(); ++x)
    out << x << ' ' << _vac_bin[x] << '\n';
}
