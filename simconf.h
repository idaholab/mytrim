#ifndef SIMCONF_H
#define SIMCONF_H

namespace MyTRIM_NS {

// ZBL coefficients a,b,c,d for all element pairs from Z=1..92
struct scoefLine {
  char sym[3], name[30];
  double mm1, m1, mnat, rho, atrho, vfermi, heat, lfctr;
};

struct simconfType {
  double ed, alfa, alpha, tmin, tau, da, cw;
  int id;

  // tables from files
  scoefLine scoef[92];
  double pcoef[92][8];
  double snuc[92][92][4];

  bool fullTraj;

  // tally electronic and nuclear stopping losses
  double EelTotal;
  double EnucTotal;

  // statistics of the simulation run
  int vacancies_created;

  simconfType( double _alfa = 0.0 );
private:
  void read_scoef();
  void read_snuc();
};

extern simconfType *simconf;

}

#endif
