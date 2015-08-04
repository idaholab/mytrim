#include "sample.h"

using namespace MyTRIM_NS;

sampleBase::sampleBase(double x, double y, double z)
{
  w[0] = x; w[1] = y; w[2] = z;
  for (unsigned int i = 0; i < 3; ++i)
    bc[i] = PBC;
}

void sampleBase::averages(const ionBase * pka)
{
  for (unsigned int i = 0; i < material.size(); ++i)
    material[i]->average(pka);
}
