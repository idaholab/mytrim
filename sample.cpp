#include "sample.h"

sampleBase::sampleBase( float x, float y, float z ) 
{ 
  w[0] = x; w[1] = y; w[2] = z;
  for( int i = 0; i < 3; i++ ) bc[i] = PBC;
}

void sampleBase::averages( const ionBase *pka )
{
  for( int i = 0; i < material.size(); i++ ) material[i]->average( pka );
}
