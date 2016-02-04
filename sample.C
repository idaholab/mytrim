#include "sample.h"

using namespace MyTRIM_NS;

SampleBase::SampleBase(Real x, Real y, Real z)
{
  w[0] = x; w[1] = y; w[2] = z;
  for (unsigned int i = 0; i < 3; ++i)
    bc[i] = PBC;
}

void
SampleBase::averages(const IonBase * pka)
{
  for (unsigned int i = 0; i < material.size(); ++i)
    material[i]->average(pka);
}

Real
SampleBase::rangeMaterial(Point & /* pos */, Point & /* dir */)
{
  return 100000.0;
}
