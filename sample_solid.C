#include "sample_solid.h"

using namespace MyTRIM_NS;

MaterialBase*  sampleSolid::lookupMaterial(Point & /* pos */)
{
  return material[0];
}
