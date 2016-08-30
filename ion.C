#include <cmath>

#include "ion.h"

using namespace MyTRIM_NS;

IonBase::IonBase() :
    _seed(0),
    _tag(-1),
    _Ef(3.0),
    _state(MOVING)
{
}

IonBase::IonBase(IonBase* prototype) :
    _Z(prototype->_Z),
    _m(prototype->_m),
    _E(prototype->_E),
    _seed(0),
    _tag(-1),
    _Ef(prototype->_Ef),
    _state(MOVING)
{
}

IonBase::IonBase(int Z, Real m, Real E) :
    _Z(Z),
    _m(m),
    _E(E),
    _seed(0),
    _tag(-1),
    _Ef(3.0),
    _state(MOVING)
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
IonBase::parent(IonBase * parent)
{
  _Ef = 3.0; // final energy

  _gen = parent->_gen + 1;
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
       << i._id << ' ' << i._gen << ' ' << i._tag << ' ';
    return os;
  }
}

bool
IonBase::operator< (const IonBase & a) const
{
  return (_Z < a._Z) || (_Z == a._Z && _m < a._m);
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

void
IonClock::parent(IonBase * parent)
{
  IonBase::parent(parent);

  IonClock * clock_parent = dynamic_cast<IonClock *>(parent);
  _time = clock_parent ? clock_parent->_time : 0.0;
}

namespace MyTRIM_NS {
  // leverage the parent class output and augment it
  std::ostream & operator << (std::ostream & os, const IonClock & i)
  {
    os << (static_cast<const IonBase &>(i)) <<  i._time << ' ';
    return os;
  }
}
