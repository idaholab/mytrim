#include "sample_dynamic.h"
#include "simconf.h"
#include "element.h"
#include <cmath>

void sampleDynamic::averages( const ionBase *_pka )
{
  // remember pka, we do not calculate averages right now, but on demand!
  pka = _pka;

  // reset update status, to trigger on-demand updates
  int i = layerUpdated.size();
  while( i >= 0 ) layerUpdated[--i] = false;
}

int sampleDynamic::lookupLayer( float* pos )
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

materialBase*  sampleDynamic::lookupMaterial( float* pos )
{
  int i = lookupLayer(pos);

  // on-demand update
  if( !layerUpdated[i] )
    material[i]->average( pka );

  return material[i];
}

void sampleDynamic::addAtomsToLayer( int layer, int n, int Z )
{
  int i, ne = material[layer]->element.size();
  double rnorm = 0.0; // volume of one atom with the fractional composition of the layer at natural density

  // look which element in the layer we are modifying and take weighted average of 1/rho of elements
  for( i = 0; i < material[layer]->element.size(); i++ )
  {
    if( material[layer]->element[i]->z == Z ) ne = i;
    rnorm += material[layer]->element[i]->t / simconf->scoef[material[layer]->element[i]->z-1].rho;
  }

  // element not yet contained in layer
  if( i == material[layer]->element.size() )
  {
    elementBase* element = new elementBase;
    element->z = Z;
    element->m = simconf->scoef[Z-1].m1;
    element->t = 0.0;
    material[layer]->element.push_back(element);
  }

  // mass of layer (g/cm^3 * Ang^3 = 10^-24 g)
  double ml = material[layer]->rho * w[1] * w[2] * layerThickness[layer];
  int nl    = material[layer]->arho * w[1] * w[2] * layerThickness[layer];

  // mass change of layer ( g/mole = g/(0.6022*10^24))
  double ma = 0.6022 * material[layer]->element[i]->m * double(n);

  // keep density consistent...
  layerThickness[layer] += ma / ( material[layer]->rho * w[1] * w[2] );

  // change stoichiometry (arho 1/ang^3 * ang^3 )
  for( i = 0; i < material[layer]->element.size(); i++ )
  {
    if( material[layer]->element[i]->z == Z )
      material[layer]->element[i]->t = ( material[layer]->element[i]->t * nl + n ) / ( nl + n );
    else
      material[layer]->element[i]->t = ( material[layer]->element[i]->t * nl ) / ( nl + n );
  }

  // mark as dirty, just in case, and re-prepare to update averages
  layerUpdated[i] = false;
  material[layer]->prepare();
}
