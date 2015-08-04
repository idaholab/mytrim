#include "sample_solid.h"

using namespace MyTRIM_NS;

materialBase*  sampleSolid::lookupMaterial(double * /* pos */)
{
  return material[0];
}
