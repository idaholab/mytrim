#ifndef THREADEDTRIMBASE_H
#define THREADEDTRIMBASE_H

#include "trim.h"

using namespace MyTRIM_NS;

class ThreadedTrimBase : public TrimBase
{
public:
  ThreadedTrimBase(SimconfType * simconf, SampleBase * sample) :
    TrimBase(simconf, sample) {}

  virtual bool followRecoil() { return !_primaries_only; }
  virtual void threadJoin(const ThreadedTrimBase &) = 0;

  bool _primaries_only;
};

#endif //THREADEDTRIMBASE_H
