#ifndef MYTRIM_ION_H
#define MYTRIM_ION_H

#include <iostream>
#include "simconf.h"

namespace MyTRIM_NS {

class IonBase
{
public:
  IonBase();
  IonBase(IonBase* prototype);
  IonBase(int Z, Real m, Real E);
  virtual ~IonBase() {};

  virtual void parent(IonBase* parent);
  virtual IonBase* spawnRecoil();

  void setEf();

  bool operator< (const IonBase &) const;

  /// atomic number
  int _Z;

  /// mass in [amu]
  Real _m;

  /// kinetic energy in [eV]
  Real _E;

  // normalized velocity vector, and position
  Point _dir, _pos;

  // random number generator seed for this ion and its recoils
  unsigned int _seed;

  // recoil generation number, unique ID
  int _gen, _id;

  // material tag
  int _tag;

  // final energy up to which this recoil will be followed in [eV]
  Real _Ef;

  // state of the recoil:
  //   MOVING          ion is still being tracked
  //   REPLACEMENT     pka->_Z == element._Z && pka->_E < element._Edisp
  //   SUBSTITUTIONAL  pka->_Z != element._Z && pka->_E < element._Edisp
  //   INTERSTITIAL    no recoil spawned and pke->_E < pka->_Ef
  //   LOST            ion has left the sample
  enum StateType { MOVING, REPLACEMENT, SUBSTITUTIONAL, INTERSTITIAL, LOST } _state;
  static const int DELETE = -1;
};

/// Serialize ion into text stream
std::ostream& operator << (std::ostream& os, const IonBase &i);


class IonMDTag : public IonBase
{
public:
  IonMDTag() : IonBase(), _md(0) {}
  IonMDTag(IonMDTag * prototype) : IonBase(prototype), _md(prototype->_md) {}

  /// overwrite this to return recoil ion objects of type IonMDTag
  virtual IonBase* spawnRecoil();

  /// generation after first ion falling into the MD energy gap (200eV - 12000eV) TODO: move to subclass?
  int _md;
};

/// Serialize ion into text stream
std::ostream& operator << (std::ostream& os, const IonMDTag &p);

class IonClock : public IonBase
{
public:
  IonClock() : IonBase(), _time(0.0) {}
  IonClock(IonClock * prototype) : IonBase(prototype), _time(prototype->_time) {}

  void parent(IonBase *parent);

  Real _time;
};

}

#endif
