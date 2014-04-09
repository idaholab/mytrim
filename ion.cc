#include <cmath>

#include "ion.h"

ionBase::ionBase()
{ 
  ef = 5.0; // final energy
  t = 0.0; //clock
}

ionBase::ionBase( ionBase* prototype )
{ 
  ef = prototype->ef; // final energy
  t = prototype->t; //clock

  z1 = prototype->z1;
  m1 = prototype->m1;
  e = prototype->e;
}

ionBase::ionBase( int _z1, double _m1, double _e ) : z1(_z1), m1(_m1), e(_e)
{
  ef = 5.0; // final energy
  t = 0.0; //clock
}

void ionBase::set_ef()
{
  ef = fmax( 5.0, 0.00001 * e ); // final energy TODO: 100Mev*0.00001 = 1keV - do we really want to stop there?!
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
