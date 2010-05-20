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
