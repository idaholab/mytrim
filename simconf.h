#ifndef SIMCONF_H
#define SIMCONF_H 1

// ZBL coefficients a,b,c,d for all element pairs from Z=1..92
struct scoefLine {
  char sym[3], name[30];
  float mm1, m1, mnat, rho, atrho, vfermi, heat, lfctr;
};

struct simconfType {
  float ed, alfa, alpha, tmin, tau, da, cw;
  int id;

  // tables from files
  scoefLine scoef[92];
  float pcoef[92][8];
  float snuc[92][92][4];

  bool fullTraj;

  simconfType( float _alfa = 0.0 );
private:
  void read_scoef();
  void read_snuc();
};

extern simconfType *simconf;

#endif