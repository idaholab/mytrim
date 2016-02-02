/*
MyTRIM - a three dimensional binary collision Monte Carlo library.
Copyright (C) 2008-2015  Daniel Schwen <daniel@schwen.de>

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA
*/

#ifndef MYTRIM_TRIM_H
#define MYTRIM_TRIM_H

#include <vector>
#include <queue>
#include <cmath>

#include "material.h"
#include "sample.h"
#include "simconf.h"

namespace MyTRIM_NS {

class trimBase {
public:
  void trim(ionBase *pka, std::queue<ionBase*> &recoils);
  trimBase(simconfType * simconf_, sampleBase *sample_) :
    simconf(simconf_), sample(sample_) {}

protected:
  simconfType * simconf;
  sampleBase *sample;
  ionBase *pka, *recoil;
  materialBase *material;
  elementBase *element;
  std::queue<ionBase*> *recoil_queue_ptr;
  bool terminate;

  // by default only follow recoils with E > 12eV
  virtual bool followRecoil() {
    // TODO: find a better place for this!
    /*if (pka->md > 0)
      recoil->md = pka->md +1;
    else
      recoil->md = 0;
    */
    return true;
  };
  virtual void vacancyCreation();
  virtual void checkPKAState() {};
  virtual void dissipateRecoilEnergy() {};
};


//
// Only follow the primary knock ons (i.e. fission fragments)
//
class trimPrimaries : public trimBase {
public:
  trimPrimaries(simconfType * simconf_, sampleBase *sample_) : trimBase(simconf_, sample_) {};
protected:
  virtual int maxGen() { return 1; };
  virtual bool followRecoil() { return (recoil->gen < maxGen()); };

  void vacancyCreation()
  {
    simconf->vacancies_created++;

    // Modified Kinchin-Pease
    if (recoil->gen == maxGen())
    {
      // calculate modified kinchin pease data
      // http://www.iue.tuwien.ac.at/phd/hoessinger/node47.html
      Real ed = 0.0115 * std::pow(material->az, -7.0/3.0) * recoil->e;
      Real g = 3.4008 * std::pow(ed, 1.0/6.0) + 0.40244 * std::pow(ed, 3.0/4.0) + ed;
      Real kd = 0.1337 * std::pow(material->az, 2.0/3.0) / std::pow(material->am, 0.5); //Z,M
      Real Ev = recoil->e / (1.0 + kd * g);
      simconf->vacancies_created += int(0.8 * Ev / (2.0*element->Edisp));

      // TODO: this is missing the energy threshold of 2.5Ed!!!!
      // TODO: should be something like material->Edisp (average?)
    }
  }
};


//
// Only follow the first generation of recoils
//
class trimRecoils : public trimPrimaries {
  public:
    trimRecoils(simconfType * simconf_, sampleBase *sample_) : trimPrimaries(simconf_, sample_) {};
  protected:
    virtual int maxGen() { return 2; };
};


//
// store a history of all recoils
//
class trimHistory : public trimBase {
public:
  trimHistory(simconfType * simconf_, sampleBase *sample_) : trimBase(simconf_, sample_) {};
  std::vector<Real> pos_hist[3];
protected:
  virtual bool followRecoil()
  {
    pos_hist[0].push_back(pka->pos(0));
    pos_hist[1].push_back(pka->pos(1));
    pos_hist[2].push_back(pka->pos(2));
    return true;
  };
};


//
// Log vaccancy/interstitial creation
//
class trimDefectLog : public trimBase {
public:
  trimDefectLog(simconfType * simconf_, sampleBase *sample_, std::ostream &os_) : trimBase(simconf_, sample_), os(os_) {};
protected:
  std::ostream &os;

  // ions being removed from lattice sites
  virtual void vacancyCreation() {
    os << "V " << *recoil << std::endl;
  };

  // ions coming to rest
  virtual void checkPKAState() {
    if (pka->state==ionBase::INTERSTITIAL)
      os << "I " << *pka << std::endl;
    else if (pka->state==ionBase::SUBSTITUTIONAL)
      os << "S " << *pka << std::endl;
    else if (pka->state==ionBase::REPLACEMENT)
      os << "R " << *pka << std::endl;
  };
};


//
// Map vaccancy creation
//
class trimVacMap : public trimBase {
  static const int mx = 20, my = 20;
public:
  int vmap[mx][my][3];
  trimVacMap(simconfType * simconf_, sampleBase *sample_, int z1_, int z2_, int z3_ = -1) : trimBase(simconf_, sample_), z1(z1_), z2(z2_), z3(z3_)
  {
    for (int e = 0; e < 3; e++)
      for (int x = 0; x < mx; x++)
        for (int y = 0; y < my; y++)
          vmap[x][y][e] = 0;
  };
protected:
  int z1, z2, z3;
  virtual void vacancyCreation()
  {
    // both atoms have enough energy to leave the site
    int x, y;

    x = ((recoil->pos(0) * mx) / sample->w[0]);
    y = ((recoil->pos(1) * my) / sample->w[1]);
    x -= int(x/mx) * mx;
    y -= int(y/my) * my;

    // keep track of vaccancies for the two constituents
    if (recoil->_Z == z1) vmap[x][y][0]++;
    else if (recoil->_Z == z2) vmap[x][y][1]++;
    else if (recoil->_Z == z3) vmap[x][y][2]++;
  };
};


//
// Output all phonon energy losses
//
class trimPhononOut : public trimBase {
public:
  trimPhononOut(simconfType * simconf_, sampleBase *sample_,  std::ostream &os_) : trimBase(simconf_, sample_), os(os_) {};
protected:
  std::ostream &os;

  // residual energy of pka coming to a stop
  virtual void checkPKAState() {
    if (pka->state == ionBase::MOVING ||
        pka->state == ionBase::LOST) return;

    os << pka->e << ' ' <<  *pka << std::endl;
    simconf->EnucTotal += pka->e;
  }

  // recoil atom is not leaving its site
  // (make sure it keeps its binding enrgy and dissipate emeining E)
  virtual void dissipateRecoilEnergy()
  {
    Real Edep = recoil->e + element->Elbind;
    os << Edep << ' ' <<  *recoil << std::endl;
    simconf->EnucTotal += Edep;
  };

  // dissipate lattice binding energy where recoil branches off
  virtual bool followRecoil() {
    os << element->Elbind << ' ' <<  *recoil << std::endl;
    simconf->EnucTotal += element->Elbind;
    return true;
  }
};

}

#endif
