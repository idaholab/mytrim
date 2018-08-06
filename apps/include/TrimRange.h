#ifndef TRIMRANGE_H
#define TRIMRANGE_H

#include "ThreadedTrimBase.h"

using namespace MyTRIM_NS;

class TrimRange : public ThreadedTrimBase
{
public:
  TrimRange(SimconfType * simconf, SampleBase * sample);

protected:
  virtual void vacancyCreation();
  virtual void dissipateRecoilEnergy();
  virtual bool followRecoil() { return (_recoil->_gen < 1); }

  virtual void threadJoin(const ThreadedTrimBase & ttb);
  virtual void writeOutput();

private:
  /// list of ranges per Z
  std::vector<std::vector<Real>> _range;
};

#endif // TRIMRANGE_H
