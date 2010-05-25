#ifndef ION_H
#define ION_H 1

struct ionBase {
  int z1;
  float m1, e;
  float dir[3], pos[3]; // normalized velocity vector, and position
  float t; // internal clock

  int tag, gen, id;
  int md; // generation after first ion falling into the MD energy gap ( 200eV - 12000eV )

  float ef;

  ionBase();
  ionBase( ionBase* parent );
  virtual ~ionBase() {};

  virtual ionBase* spawnRecoil();

  void set_ef();
};

#endif
