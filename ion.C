#include <cmath>

#include "ion.h"

using namespace MyTRIM_NS;

IonBase::IonBase() :
    _time(0.0),
    _Ef(3.0),
    state(MOVING)
{
}

IonBase::IonBase(IonBase* prototype) :
    _Z(prototype->_Z),
    _m(prototype->_m),
    _E(prototype->_E),
    _time(prototype->_time),
    _Ef(prototype->_Ef),
    state(MOVING)
{
}

IonBase::IonBase(int Z, Real m, Real E) :
    _Z(Z),
    _m(m),
    _E(E),
    _time(0.0),
    _Ef(3.0),
    state(MOVING)
{
}

void
IonBase::setEf()
{
  // stop following an ion if it's energy falls below 5.0eV
  _Ef = 3.0;

  // final energy TODO: 100Mev*0.00001 = 1keV - do we really want to stop there?!
  //fmax(5.0, 0.00001 * e);
}

void
IonBase::parent(IonBase *parent)
{
  _Ef = 3.0; // final energy

  gen = parent->gen + 1;
  _time = parent->_time;

  _pos = parent->_pos;
}

IonBase*
IonBase::spawnRecoil()
{
  IonBase * recoil = new IonBase;
  recoil->parent(this);
  return recoil;
}

// output operator (implement for derived classes if necessary)
namespace MyTRIM_NS {
  std::ostream & operator << (std::ostream & os, const IonBase & i)
  {
    os << i._pos(0) << ' ' << i._pos(1) << ' ' << i._pos(2) << ' '
       << i._Z << ' ' << i._m << ' ' << i._E << ' '
       << i._time << ' '
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
  std::ostream & operator << (std::ostream & os, const IonMDTag & i)
  {
    os << (static_cast<const IonBase &>(i)) <<  i._md << ' ';
    return os;
  }
}
