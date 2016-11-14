#ifndef THREADEDTRIMBASE_H
#define THREADEDTRIMBASE_H

#include "trim.h"

using namespace MyTRIM_NS;

class ThreadedTrimBase : public TrimBase
{
public:
  ThreadedTrimBase(SimconfType * simconf, SampleBase * sample) :
    TrimBase(simconf, sample) {}

  virtual void threadJoin(const ThreadedTrimBase &) = 0;
};

#endif //THREADEDTRIMBASE_H
