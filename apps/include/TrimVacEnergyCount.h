#ifndef TRIMVACENERGYCOUNT_H
#define TRIMVACENERGYCOUNT_H

#include "ThreadedTrimBase.h"

using namespace MyTRIM_NS;

class TrimVacEnergyCount : public ThreadedTrimBase
{
public:
  TrimVacEnergyCount(SimconfType * simconf, SampleBase * sample);

protected:
  virtual void vacancyCreation();

  virtual void threadJoin(const ThreadedTrimBase & ttb);
  virtual void writeOutput();

private:
  /// histogram of vacancies created per log10 Energy and unit depth (Ang)
  std::vector<std::vector<unsigned int> > _evac_bin;
};

#endif //TRIMVACENERGYCOUNT_H
