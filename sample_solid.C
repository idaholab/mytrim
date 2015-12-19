#include "sample_solid.h"

using namespace MyTRIM_NS;

materialBase*  sampleSolid::lookupMaterial(Real * /* pos */)
{
  return material[0];
}
