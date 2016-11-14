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

#ifndef MYTRIM_SIMCONF_H
#define MYTRIM_SIMCONF_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <memory>

#ifdef MYTRIM_ENABLED
// building from within MOOSE/Magpie
#include "MooseError.h"
#include "MooseTypes.h"
#include "libmesh/point.h"
#include "libmesh/utility.h"
#else
// building standalone (for Travis CI tests)
typedef double Real;
#include "shim/pow.h"
#include "shim/point.h"
#endif

namespace MyTRIM_NS {

class SimconfType
{
public:
  SimconfType(unsigned int seed = 12345678);

  inline Real drand() { return cxx11random_dis_Real(*cxx11random_gen); }
  inline unsigned int irand() { return cxx11random_dis_int(*cxx11random_gen); }
  void seed(unsigned int seed);

  // length scale in Angstrom
  // 1  =>    Units are Angstrom
  // 10 =>    Units are Nanometers
  // 10000 => Units are Mircometers
  // etc.
  void setLengthScale(Real l);
  Real lengthScale() { return _length_scale; }
  Real areaScale() { return _area_scale; }
  Real volumeScale() { return _volume_scale; }

  Real ed, tmin, tau, da, cw;
  int _id;

  // tables from files
  // ZBL coefficients a, b, c, d for all element pairs from Z=1..92
  struct ScoefLine
  {
    ScoefLine();
    void read95A(std::ifstream &);
    void read95B(std::ifstream &);
    void readSlfctr(std::ifstream &);
    void readElname(std::ifstream &);

    // from ELNAME.dat
    std::string sym, name;

    // from SCOEF.95A
    Real mm1, m1, mnat, rho, atrho, vfermi, heat, lfctr;
    std::vector<Real> pcoef;     // 8 columns

    // from SCOEF.95B
    std::vector<Real> ehigh;     // 4 columns
    std::vector<Real> screen;    // 19 columns
    std::vector<Real> fermicorr; // 15 columns
  };

  // tables from files
  static const unsigned int _rows = 92;
  std::vector<ScoefLine> scoef;
  ScoefLine scoeflast;

  /// ZBL coefficients
  Real snuc[_rows][_rows][4];

  /// dump data for full trajectories to std::cout
  bool fullTraj;

  ///@{ tally electronic and nuclear stopping losses
  Real EelTotal;
  Real EnucTotal;
  ///@}

  /// statistics of the simulation run
  int vacancies_created;

private:
  void readDataFiles();

  static void fileReadError(const std::string &);

  void skipLine(FILE * sf);
  void skipLine(std::ifstream & sf);

  std::string _data_dir;

  std::unique_ptr<std::mt19937> cxx11random_gen;
  std::uniform_real_distribution<double> cxx11random_dis_Real;
  std::uniform_int_distribution<unsigned int> cxx11random_dis_int;

  Real _length_scale;
  Real _area_scale;
  Real _volume_scale;
};

}

#endif
