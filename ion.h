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
  IonBase(int Z, Real m, Real e_);
  virtual ~IonBase() {};

  virtual void parent(IonBase* parent);
  virtual IonBase* spawnRecoil();

  void setEf();

  /// atomic number
  int _Z;

  /// mass
  Real _m;

  /// kinetic energy
  Real e;

  // normalized velocity vector, and position
  Point dir, pos;

  // internal clock (needed for visualization)
  Real t;

  // recoil generation number, unique ID
  int gen, id;

  // material tag
  int tag;

  // final energy up to which this recoil will be followed
  Real ef;

  // state of the recoil:
  //   MOVING          ion is still being tracked
  //   REPLACEMENT     pka->_Z == element->_Z && pka->e < element->_Edisp
  //   SUBSTITUTIONAL  pka->_Z != element->_Z && pka->e < element->_Edisp
  //   INTERSTITIAL    no recoil spawned and pke->e < pka->ef
  //   LOST            ion has left the sample
  enum StateType { MOVING, REPLACEMENT, SUBSTITUTIONAL, INTERSTITIAL, LOST } state;
  static const int DELETE = -1;
};

/// Serialize ion into text stream
std::ostream& operator << (std::ostream& os, const IonBase &i);


class IonMDTag : public IonBase
{
public:
  IonMDTag() : IonBase(), md(0) {}
  IonMDTag(IonMDTag * prototype) : IonBase(prototype), md(prototype->md) {}

  /// overwrite this to return recoil ion objects of type IonMDTag
  virtual IonBase* spawnRecoil();

  /// generation after first ion falling into the MD energy gap (200eV - 12000eV) TODO: move to subclass?
  int md;
};

/// Serialize ion into text stream
std::ostream& operator << (std::ostream& os, const IonMDTag &p);
}

#endif
