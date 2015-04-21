#ifndef ION_H
#define ION_H

#include <iostream>

namespace MyTRIM_NS {

struct ionBase {
  // atomic number, mass, and kinetic energy of the ion
  int z1;
  double m1, e;

  // normalized velocity vector, and position
  double dir[3], pos[3];

  // internal clock (needed for visualization)
  double t;

  // recoil generation number, unique ID
  int gen, id;

  // material tag
  int tag;

  // final energy up to which this recoil will be followed
  double ef;

  // state of the recoil:
  //   MOVING          ion is still being tracked
  //   REPLACEMENT     pka->z1 == element->z && pka->e < element->Edisp
  //   SUBSTITUTIONAL  pka->z1 != element->z && pka->e < element->Edisp
  //   INTERSTITIAL    no recoil spawned and pke->e < pka->ef
  //   LOST            ion has left the sample
  enum StateType { MOVING, REPLACEMENT, SUBSTITUTIONAL, INTERSTITIAL, LOST } state;
  static const int DELETE = -1;

  ionBase();
  ionBase( ionBase* prototype );
  ionBase( int _z1, double _m1, double _e );
  virtual ~ionBase() {};

  virtual void parent( ionBase* parent );
  virtual ionBase* spawnRecoil();

  void set_ef();
};

std::ostream& operator << (std::ostream& os, const ionBase &i); /// Serialize ion into text stream


struct ionMDtag : public ionBase {
  // generation after first ion falling into the MD energy gap ( 200eV - 12000eV ) TODO: move to subclass?
  int md;

  // overwrite this tor return recoil ion objects of type ionMDtag
  virtual ionBase* spawnRecoil();
};

std::ostream& operator << (std::ostream& os, const ionMDtag &p); /// Serialize ion into text stream

}

#endif
