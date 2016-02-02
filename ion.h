#ifndef MYTRIM_ION_H
#define MYTRIM_ION_H

#include <iostream>
#include "simconf.h"

namespace MyTRIM_NS {

class ionBase
{
public:
  ionBase();
  ionBase(ionBase* prototype);
  ionBase(int Z, Real m, Real e_);
  virtual ~ionBase() {};

  virtual void parent(ionBase* parent);
  virtual ionBase* spawnRecoil();

  void set_ef();

  /// atomic number
  int _Z;

  /// mass
  Real _m;

  /// kinetic energy
  Real e;

  // normalized velocity vector, and position
  Real dir[3], pos[3];

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
  //   REPLACEMENT     pka->_Z == element->z && pka->e < element->Edisp
  //   SUBSTITUTIONAL  pka->_Z != element->z && pka->e < element->Edisp
  //   INTERSTITIAL    no recoil spawned and pke->e < pka->ef
  //   LOST            ion has left the sample
  enum StateType { MOVING, REPLACEMENT, SUBSTITUTIONAL, INTERSTITIAL, LOST } state;
  static const int DELETE = -1;
};

/// Serialize ion into text stream
std::ostream& operator << (std::ostream& os, const ionBase &i);


class ionMDtag : public ionBase
{
public:
  ionMDtag() : ionBase(), md(0) {}
  ionMDtag(ionMDtag * prototype) : ionBase(prototype), md(prototype->md) {}

  /// overwrite this to return recoil ion objects of type ionMDtag
  virtual ionBase* spawnRecoil();

  /// generation after first ion falling into the MD energy gap (200eV - 12000eV) TODO: move to subclass?
  int md;
};

/// Serialize ion into text stream
std::ostream& operator << (std::ostream& os, const ionMDtag &p);
}

#endif
