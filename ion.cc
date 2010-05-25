#include <math.h>

#include "ion.h"

ionBase::ionBase()
{ 
  ef = 5.0; // final energy
  t = 0.0; //clock
}

void ionBase::set_ef()
{ 
  ef = fmax( 5.0, 0.00001 * e ); // final energy
}

void ionBase::parent( ionBase *parent )
{
  ef = 5.0; // final energy

  gen = parent->gen + 1;
  t = parent->t;

  for( int i = 0; i < 3; i++ ) 
    pos[i] = parent->pos[i];
}

ionBase* ionBase::spawnRecoil()
{
  ionBase *recoil = new ionBase;
  recoil->parent(this);
  return recoil;
}
