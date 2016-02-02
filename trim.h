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

class TrimBase
{
public:
  TrimBase(SimconfType * simconf_, SampleBase *sample_) :
      simconf(simconf_),
      sample(sample_)
  {}

  void trim(IonBase *pka, std::queue<IonBase*> &recoils);

protected:
  SimconfType * simconf;
  SampleBase * sample;
  IonBase * pka, * recoil;
  MaterialBase * material;
  ElementBase * element;
  std::queue<IonBase*> * recoil_queue_ptr;
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
  virtual void checkPKAState() {}
  virtual void dissipateRecoilEnergy() {}
};


//
// Only follow the primary knock ons (i.e. fission fragments)
//
class TrimPrimaries : public TrimBase
{
public:
  TrimPrimaries(SimconfType * simconf, SampleBase * sample) :
      TrimBase(simconf, sample)
  {}

protected:
  virtual int maxGen() { return 1; }
  virtual bool followRecoil() { return (recoil->gen < maxGen()); }

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
      simconf->vacancies_created += int(0.8 * Ev / (2.0*element->_Edisp));

      // TODO: this is missing the energy threshold of 2.5Ed!!!!
      // TODO: should be something like material->_Edisp (average?)
    }
  }
};


//
// Only follow the first generation of recoils
//
class TrimRecoils : public TrimPrimaries
{
public:
  TrimRecoils(SimconfType * simconf, SampleBase * sample) :
      TrimPrimaries(simconf, sample)
  {}

protected:
  virtual int maxGen() { return 2; }
};


//
// store a history of all recoils
//
class TrimHistory : public TrimBase
{
public:
  TrimHistory(SimconfType * simconf, SampleBase * sample) :
      TrimBase(simconf, sample)
  {}

  const std::vector<Point> & getHistory() { return _pos_hist; }

protected:
  virtual bool followRecoil()
  {
    _pos_hist.push_back(pka->pos);
    return true;
  }

  std::vector<Point> _pos_hist;
};


//
// Log vaccancy/interstitial creation
//
class TrimDefectLog : public TrimBase
{
public:
  TrimDefectLog(SimconfType * simconf, SampleBase * sample, std::ostream & os) :
      TrimBase(simconf, sample),
      _os(os)
  {}

protected:
  std::ostream & _os;

  /// ions being removed from lattice sites
  virtual void vacancyCreation() {
    _os << "V " << *recoil << std::endl;
  }

  /// ions coming to rest
  virtual void checkPKAState() {
    if (pka->state==IonBase::INTERSTITIAL)
      _os << "I " << *pka << std::endl;
    else if (pka->state==IonBase::SUBSTITUTIONAL)
      _os << "S " << *pka << std::endl;
    else if (pka->state==IonBase::REPLACEMENT)
      _os << "R " << *pka << std::endl;
  }
};


//
// Map vaccancy creation
//
class TrimVacMap : public TrimBase
{
  static const int mx = 20, my = 20;

public:
  TrimVacMap(SimconfType * simconf, SampleBase * sample, int z1, int z2, int z3 = -1) :
      TrimBase(simconf, sample),
      _z1(z1),
      _z2(z2),
      _z3(z3)
  {
    for (int e = 0; e < 3; e++)
      for (int x = 0; x < mx; x++)
        for (int y = 0; y < my; y++)
          vmap[x][y][e] = 0;
  }

  int vmap[mx][my][3];

protected:
  // TODO: this is horrible. Pass in variable number of Z!
  int _z1, _z2, _z3;

  virtual void vacancyCreation()
  {
    // both atoms have enough energy to leave the site
    int x, y;

    x = ((recoil->pos(0) * mx) / sample->w[0]);
    y = ((recoil->pos(1) * my) / sample->w[1]);
    x -= int(x/mx) * mx;
    y -= int(y/my) * my;

    // keep track of vaccancies for the two constituents
    if (recoil->_Z == _z1)
      vmap[x][y][0]++;
    else if (recoil->_Z == _z2)
      vmap[x][y][1]++;
    else if (recoil->_Z == _z3)
      vmap[x][y][2]++;
  }
};


//
// Output all phonon energy losses
//
class TrimPhononOut : public TrimBase
{
public:
  TrimPhononOut(SimconfType * simconf, SampleBase * sample,  std::ostream & os) :
      TrimBase(simconf, sample),
      _os(os)
  {}

protected:
  std::ostream & _os;

  // residual energy of pka coming to a stop
  virtual void checkPKAState()
  {
    if (pka->state == IonBase::MOVING ||
        pka->state == IonBase::LOST) return;

    _os << pka->e << ' ' <<  *pka << std::endl;
    simconf->EnucTotal += pka->e;
  }

  // recoil atom is not leaving its site
  // (make sure it keeps its binding enrgy and dissipate emeining E)
  virtual void dissipateRecoilEnergy()
  {
    Real Edep = recoil->e + element->_Elbind;
    _os << Edep << ' ' <<  *recoil << std::endl;
    simconf->EnucTotal += Edep;
  }

  // dissipate lattice binding energy where recoil branches off
  virtual bool followRecoil() {
    _os << element->_Elbind << ' ' <<  *recoil << std::endl;
    simconf->EnucTotal += element->_Elbind;
    return true;
  }
};

}

#endif
