#include "sample_solid.h"

using namespace MyTRIM_NS;

materialBase*  sampleSolid::lookupMaterial(Point & /* pos */)
{
  return material[0];
}
