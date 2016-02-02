#include <cmath>

#include "ion.h"

using namespace MyTRIM_NS;

IonBase::IonBase() :
    t(0.0),  // clock
    ef(3.0), // final energy
    state(MOVING)
{
}

IonBase::IonBase(IonBase* prototype) :
    _Z(prototype->_Z),
    _m(prototype->_m),
    e(prototype->e),
    state(MOVING)
{
  ef = prototype->ef; // final energy
  t = prototype->t;   //clock
}

IonBase::IonBase(int Z, Real m, Real e_) :
    _Z(Z),
    _m(m),
    e(e_),
    state(MOVING)
{
  ef = 3.0; // final energy
  t = 0.0; //clock;
}

void
IonBase::setEf()
{
  // stop following an ion if it's energy falls below 5.0eV
  ef = 3.0;

  // final energy TODO: 100Mev*0.00001 = 1keV - do we really want to stop there?!
  //fmax(5.0, 0.00001 * e);
}

void
IonBase::parent(IonBase *parent)
{
  ef = 3.0; // final energy

  gen = parent->gen + 1;
  t = parent->t;

  for (unsigned int i = 0; i < 3; i++)
    pos(i) = parent->pos(i);
}

IonBase*
IonBase::spawnRecoil()
{
  IonBase *recoil = new IonBase;
  recoil->parent(this);
  return recoil;
}

// output operator (implement for derived classes if necessary)
namespace MyTRIM_NS {
  std::ostream& operator << (std::ostream& os, const IonBase &i)
  {
    os << i.pos(0) << ' ' << i.pos(1) << ' ' << i.pos(2) << ' '
       << i._Z << ' ' << i._m << ' ' << i.e << ' '
       << i.t << ' '
       << i.id << ' ' << i.gen << ' ' << i.tag << ' ';
    return os;
  }
}


IonBase*
IonMDTag::spawnRecoil()
{
  IonBase *recoil = new IonMDTag;
  recoil->parent(this);
  return recoil;
}

namespace MyTRIM_NS {
  // leverage the parent class output and augment it
  std::ostream& operator << (std::ostream& os, const IonMDTag &i)
  {
    os << (static_cast<const IonBase &>(i)) <<  i.md << ' ';
    return os;
  }
}
