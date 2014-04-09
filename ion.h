#ifndef ION_H
#define ION_H 1

struct ionBase {
  // atomic number, mass, and kinetic energy of the ion
  int z1;
  double m1, e;

  // normalized velocity vector, and position
  double dir[3], pos[3]; 

  // internal clock (needed for visualization)
  double t; 

  // integer tag, recoild generation number, ID field
  int tag, gen, id;

  // generation after first ion falling into the MD energy gap ( 200eV - 12000eV ) TODO: move to subclass?
  int md; 

  // final energy up to which this recoil will be followed
  double ef;

  // state of the recoil
  enum StateType { MOVING, REPLACEMENT, INTERSTITIAL } state;

  ionBase();
  ionBase( ionBase* prototype );
  ionBase( int _z1, double _m1, double _e );
  virtual ~ionBase() {};

  virtual void parent( ionBase* parent );
  virtual ionBase* spawnRecoil();

  void set_ef();
};

#endif
