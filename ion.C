#include <cmath>

#include "ion.h"

using namespace MyTRIM_NS;

ionBase::ionBase() :
    t(0.0),  // clock
    ef(3.0), // final energy
    state(MOVING)
{
}

ionBase::ionBase(ionBase* prototype) :
    _Z(prototype->_Z),
    _m(prototype->_m),
    e(prototype->e),
    state(MOVING)
{
  ef = prototype->ef; // final energy
  t = prototype->t;   //clock
}

ionBase::ionBase(int Z, Real m, Real e_) :
    _Z(Z),
    _m(m),
    e(e_),
    state(MOVING)
{
  ef = 3.0; // final energy
  t = 0.0; //clock;
}

void ionBase::set_ef()
{
  // stop following an ion if it's energy falls below 5.0eV
  ef = 3.0;

  // final energy TODO: 100Mev*0.00001 = 1keV - do we really want to stop there?!
  //fmax(5.0, 0.00001 * e);
}

void ionBase::parent(ionBase *parent)
{
  ef = 3.0; // final energy

  gen = parent->gen + 1;
  t = parent->t;

  for (unsigned int i = 0; i < 3; i++)
    pos(i) = parent->pos(i);
}

ionBase* ionBase::spawnRecoil()
{
  ionBase *recoil = new ionBase;
  recoil->parent(this);
  return recoil;
}

// output operator (implement for derived classes if necessary)
namespace MyTRIM_NS {
  std::ostream& operator << (std::ostream& os, const ionBase &i)
  {
    os << i.pos(0) << ' ' << i.pos(1) << ' ' << i.pos(2) << ' '
       << i._Z << ' ' << i._m << ' ' << i.e << ' '
       << i.t << ' '
       << i.id << ' ' << i.gen << ' ' << i.tag << ' ';
    return os;
  }
}


ionBase* ionMDtag::spawnRecoil()
{
  ionBase *recoil = new ionMDtag;
  recoil->parent(this);
  return recoil;
}

namespace MyTRIM_NS {
  // leverage the parent class output and augment it
  std::ostream& operator << (std::ostream& os, const ionMDtag &i)
  {
    os << (static_cast<const ionBase &>(i)) <<  i.md << ' ';
    return os;
  }
}
