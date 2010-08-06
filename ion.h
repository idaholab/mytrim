#ifndef ION_H
#define ION_H 1

struct ionBase {
  int z1;
  double m1, e;
  double dir[3], pos[3]; // normalized velocity vector, and position
  double t; // internal clock

  int tag, gen, id;
  int md; // generation after first ion falling into the MD energy gap ( 200eV - 12000eV )

  double ef;

  ionBase();
  virtual ~ionBase() {};

  virtual void parent( ionBase* parent );
  virtual ionBase* spawnRecoil();

  void set_ef();
};

#endif
