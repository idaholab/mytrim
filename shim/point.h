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

#ifndef MYTRIM_POINT_H
#define MYTRIM_POINT_H

#include "simconf.h"
#include <cmath>

class Point
{
public:
  Point();
  Point(Real x, Real y, Real z);

  /// component access for backwards compatibility
  Real & operator() (unsigned int i);
  const Real & operator() (unsigned int i) const;

  /// distance squared form the origin
  Real size_sq() { return data[0]*data[0] + data[1]*data[1] + data[2]*data[2]; }

  /// distance form the origin
  Real size() { return std::sqrt(size_sq()); }

  /// arithmetic operators
  Point operator+ (const Point & rhs);
  Point operator- (const Point & rhs);
  Point operator* (Real rhs);
  Point operator/ (Real rhs);

  /// compound operators
  Point & operator+= (const Point & rhs);
  Point & operator-= (const Point & rhs);
  Point & operator*= (Real rhs);
  Point & operator/= (Real rhs);

private:
  Real data[3];
};

#endif //MYTRIM_POINT_H
