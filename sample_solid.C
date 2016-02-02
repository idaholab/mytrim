#include "sample_solid.h"

using namespace MyTRIM_NS;

MaterialBase*
SampleSolid::lookupMaterial(Point & /* pos */)
{
  return material[0];
}
