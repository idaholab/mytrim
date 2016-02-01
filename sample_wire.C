#include <cmath>
#if !defined(__APPLE__)
#include "malloc.h"
#endif

#include "sample_wire.h"

using namespace MyTRIM_NS;

sampleWire::sampleWire(Real x, Real y, Real z)  : sampleBase(x, y, z)
{
 bc[0] = CUT;
 bc[1] = CUT;
}

// look if we are within dr of the wire axis
materialBase* sampleWire::lookupMaterial(Point & pos)
{
  Real x = (pos[0] / w[0]) * 2.0 - 1.0;
  Real y = (pos[1] / w[1]) * 2.0 - 1.0;

  if ((x*x + y*y) > 1.0) return 0;
  return material[0];
}
