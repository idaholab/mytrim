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
  // cover layer
  if( pos[2] < 0.0 && pos[2] >= 250.0 )
    return material[1];

  // above sample or inside substrate
  if( pos[2] > w[2] || pos[2] < 250.0  )
    return 0;

  // in wire layer
  materialBase *ret = sampleWire::lookupMaterial(pos);
  if( ret == 0 )
    return material[1];
  else
    return ret;
}
