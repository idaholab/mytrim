#include <cmath>
#include "malloc.h"

#include "sample_wire.h"


sampleWire::sampleWire( double x, double y, double z )  : sampleBase( x, y, z)
{
 bc[0] = CUT;
 bc[1] = CUT;
}

// look if we are within dr of the wire axis
materialBase* sampleWire::lookupMaterial( double* pos ) 
{
  double x = ( pos[0] / w[0] ) * 2.0 - 1.0;
  double y = ( pos[1] / w[1] ) * 2.0 - 1.0;

  if( ( x*x + y*y ) > 1.0 ) return 0;
  return material[0];
}
