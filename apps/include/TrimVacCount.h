#ifndef TRIMVACCOUNT_H
#define TRIMVACCOUNT_H

#include "trim.h"

using namespace MyTRIM_NS;

class TrimVacCount : public TrimBase
{
public:
  TrimVacCount(SimconfType * simconf, SampleBase * sample);

protected:
  virtual void vacancyCreation();

  virtual void startOutput();
  virtual void stopOutput();

private:
  /// histogram of vacancies created per unit depth (Ang)
  std::vector<unsigned int> _vac_bin;
};

#endif //TRIMVACCOUNT_H
