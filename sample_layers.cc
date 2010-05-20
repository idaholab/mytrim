#include "sample_layers.h"

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
    return material[ material.size() - 1 ];
  else
    return material[i];
}

float sampleLayers::rangeMaterial( float* pos, float* dir )
{
  // assume dir is a normalized vector
  float range;

  if( dir[0] != 0.0 ) 
  {
    range = ( 100.0 - pos[0] ) / dir[0];
    if( range > 0 )
      return range;
  }

  return 100000.0;
}