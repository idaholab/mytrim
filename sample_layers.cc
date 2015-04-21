#include "sample_layers.h"
#include <cmath>

using namespace MyTRIM_NS;

int sampleLayers::lookupLayer( double* pos )
{
  int i;
  double d = 0.0;

  for( i = 0; i < layerThickness.size(); i++ )
  {
    d += layerThickness[i];
    if( pos[0] < d ) break;
  }

  if( i >= material.size() )
    i = material.size() - 1 ; // or 0, but we leave that to be determined by bc[] == CUT

  return i;
}

materialBase*  sampleLayers::lookupMaterial( double* pos )
{
  return material[lookupLayer(pos)];
}


double sampleLayers::rangeMaterial( double* pos, double* dir )
{
  // assume dir is a normalized vector
  double range, d = 0.0;
  int i;
  const double unrestricted = 1.0e6;

  // parallel to layer interfaces
  if( dir[0] == 0.0 ) return unrestricted;

  double epsilon = std::abs( 1.0e-10 / dir[0] );

  // outside film
  if( pos[0] < 0.0 )
  {
    if( dir[0] < 0.0 )
      return unrestricted;
    else
      return -pos[0]/dir[0] + epsilon;
  }

  // find layer
  for( i = 0; i < layerThickness.size(); i++ )
  {
    if( pos[0] >= d && pos[0] < d + layerThickness[i] )
    {
      if( dir[0] < 0 )
        return (d-pos[0])/dir[0] + epsilon;
      else
        return (d+layerThickness[i]-pos[0])/dir[0] + epsilon;
    }
    d += layerThickness[i];
  }

  // not returned yet, means we are beyond the last layer
  if( dir[0] > 0.0 )
    return unrestricted;
  else
    return (d-pos[0])/dir[0] + epsilon;
}
