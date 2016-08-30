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
#include <string>

#include "material.h"
#include "sample.h"
#include "simconf.h"

namespace MyTRIM_NS {

class TrimBase
{
public:
  TrimBase(SimconfType * simconf, SampleBase * sample) :
      _potential(UNIVERSAL),
      _simconf(simconf),
      _sample(sample),
      _base_name("mytrim"),
      _outputting(false)
  {}

  /**
   * The virtual destructor should handle closing output files
   */
  virtual ~TrimBase() { stopOutput(); }

  /**
   * Run a TRIM simulation with a given PKA and push the resulting recoils onto
   * the recoils queue
   */
  void trim(IonBase * _pka, std::queue<IonBase*> & recoils);

  /**
   * Set the output file base name
   */
  void setBaseName(const std::string & name) { _base_name = name; }

  /// overload and call baseclass version from here. Open files necessary for output in this method.
  virtual void startOutput() { _outputting = true; }

  /// overload and call baseclass version from here. Close files necessary for output in this method.
  virtual void stopOutput() { _outputting = false; }

  /// Scattering potential type
  enum Potential { UNIVERSAL, MOLIERE, CKR };
  Potential _potential;

protected:
  /// by default only follow recoils with E > 12eV
  virtual bool followRecoil();

  /// called whenever a vaccancy is created
  virtual void vacancyCreation();

  virtual void checkPKAState() {}

  /// called if recoil energy needs to get dissipated, to record phonons
  virtual void dissipateRecoilEnergy() {}

  /// helper function to determine if the output has been started
  bool outputting() { return _outputting; }

  SimconfType * _simconf;
  SampleBase * _sample;

  /// the current PKA and the last recoil it created
  IonBase * _pka, * _recoil;
  MaterialBase * _material;
  ElementBase * _element;
  std::queue<IonBase*> * recoil_queue_ptr;
  bool terminate;

  /// current path segment length
  Real _ls;

  /// current electronic energy loss along _ls
  Real _dee;

  /// current electronic energy loss at the collison after _ls
  Real _den;

  /// TRIM classes that output stuff use this string as the base name
  std::string _base_name;

private:
  /// has the output been initialized
  bool _outputting;
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
  virtual bool followRecoil() { return (_recoil->_gen < maxGen()); }
  virtual void vacancyCreation();
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
    _pos_hist.push_back(_pka->_pos);
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
  virtual void vacancyCreation();

  /// log ions coming to rest
  virtual void checkPKAState();
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
    for (unsigned int e = 0; e < 3; ++e)
      for (unsigned int x = 0; x < mx; ++x)
        for (unsigned int y = 0; y < my; ++y)
          vmap[x][y][e] = 0;
  }

  int vmap[mx][my][3];

protected:
  // TODO: this is horrible. Pass in variable number of Z!
  int _z1, _z2, _z3;

  virtual void vacancyCreation();
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

  // residual energy of _pka coming to a stop
  virtual void checkPKAState();

  // recoil atom is not leaving its site
  // (make sure it keeps its binding enrgy and dissipate emeining E)
  virtual void dissipateRecoilEnergy();

  // dissipate lattice binding energy where recoil branches off
  virtual bool followRecoil();
};

}

#endif
