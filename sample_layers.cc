#include "sample_layers.h"
#include <cmath>

materialBase*  sampleLayers::lookupMaterial( float* pos )
{
  int i;
  double d = 0.0;

  for( i = 0; i < layerThickness.size(); i++ )
  {
    d += layerThickness[i];
    if( pos[0] < d ) break;
  }

  if( i >= material.size() )
    return material[ material.size() - 1 ]; // or 0, but we leave that to be determined by bc[] == CUT
  else
    return material[i];
}

float sampleLayers::rangeMaterial( float* pos, float* dir )
{
  // assume dir is a normalized vector
  float range;

  if( dir[0] != 0.0 ) 
  {
    range = ( 500.0 - pos[0] ) / dir[0];
    if( range > 0.0 )
    {
      return range + fabs( 1.0e-10 / dir[0] ); // make it end up _just_ in the next material
    }
  }

  return 100000.0;
}