#include <cmath>
#include "malloc.h"

#include "sample_burried_wire.h"


sampleBurriedWire::sampleBurriedWire( double x, double y, double z )  : sampleWire( x, y, z)
{
 bc[0] = INF;
 bc[1] = INF;
 bc[2] = INF;
}

// look if we are within dr of the wire axis
materialBase* sampleBurriedWire::lookupMaterial( double* pos ) 
{
  materialBase *ret = sampleWire::lookupMaterial(pos);
  if( ret == 0 )
  {
    if( pos[0] > -250.0 )
      return material[1];
    else
      return 0;
  }
  else
    return ret;
}
